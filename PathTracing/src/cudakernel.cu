#include <glm/gtx/constants.hpp>

#include <curand_kernel.h>

#include "cudakernel.cuh"

__constant__ curandState_t* state;

__constant__ float PI;
__constant__ float EPSILON = 0.001f;

__constant__ int resX, resY;
__constant__ int maxDepth;
__constant__ int samples;

__constant__ float _camPos[3], _camDir[3], _camUp[3];
__constant__ float camFocal, camFovy;

__device__ GPUBVHNode* bvh = 0;
__constant__ int bvhSize = 0;

struct GPUImage
{
	int width;
	int height;
	float* data;
	__host__ __device__ GPUImage() :
		width(0),
		height(0),
		data(0)
	{
	}
};

__device__ GPUImage* textures = 0;
int numTextures = 0;

struct GPUPolarInfo
{
	glm::vec3 normal;
	glm::vec3 indir;
	glm::vec3 outdir;

	glm::vec2 R;
	glm::vec3 S;

	__device__ GPUPolarInfo()
	{
		normal = glm::vec3();
		indir = glm::vec3();
		outdir = glm::vec3();

		R = glm::vec2();
		S = glm::vec3();
	}
};

__device__ glm::vec2 GetRsRp(float n1, float n2, float cos_i)
{
	cos_i = fabs(cos_i);
	float cos_t = sqrtf(1.0f - (n1 * n1) / (n2 * n2) * (1.0f - cos_i * cos_i));
	float rs = (n1 * cos_i - n2 * cos_t) / (n1 * cos_i + n2 * cos_t);
	float rp = (n2 * cos_i - n1 * cos_t) / (n1 * cos_t + n2 * cos_i);
	return glm::vec2(rs, rp);
}

__device__ glm::mat3 GetRotationMatrix(float phi)
{
	glm::mat3 M;
	float sin_2phi = sinf(2.0f * phi);
	float cos_2phi = cosf(2.0f * phi);
	M[0] = glm::vec3(1.0f, 0.0f, 0.0f);
	M[1] = glm::vec3(0.0f, cos_2phi, -sin_2phi);
	M[2] = glm::vec3(0.0f, sin_2phi, cos_2phi);
	return M;
}

__device__ glm::vec3 CalculatePolarResult
(
	const glm::vec3& initDir,
	GPUPolarInfo* polarInfoList,
	int listSize,
	curandState_t& state
)
{
	int i = listSize;
	if (i == 0)
		return glm::vec3(0.0f);

	i--;
	glm::vec3 res = glm::vec3(0.0f);
	glm::vec3 s = glm::normalize(glm::vec3
	(
		curand_uniform(&state) * 2.0f - 1.0f,
		curand_uniform(&state) * 2.0f - 1.0f,
		curand_uniform(&state) * 2.0f - 1.0f
	));
	float phi = 2.0f * PI * curand_uniform(&state);
	for (; i >= 0; i--)
	{
		glm::vec3 d = polarInfoList[i].outdir;
		glm::vec3 s_old = s;
		s = glm::normalize(glm::cross(polarInfoList[i].normal, d));
		//phi = acos(glm::dot(s, s_old));
		phi = acos(glm::dot(s, s_old)) * glm::sign(glm::dot(glm::cross(s_old, s), d));

		res = GetRotationMatrix(phi) * res;
		const float Rs = polarInfoList[i].R.x * polarInfoList[i].R.x;
		const float Rp = polarInfoList[i].R.y * polarInfoList[i].R.y;
		glm::mat3 Refl = glm::mat3
		(
			glm::vec3((Rs + Rp) / 2.0f, (Rs - Rp) / 2.0f, 0.0f),
			glm::vec3((Rs - Rp) / 2.0f, (Rs + Rp) / 2.0f, 0.0f),
			glm::vec3(0.0f, 0.0f, polarInfoList[i].R.x * polarInfoList[i].R.y)
		);
		glm::mat3 Trans = glm::mat3
		(
			glm::vec3(1.0f - (Rs + Rp) / 2.0f, (Rp - Rs) / 2.0f, 0.0f),
			glm::vec3((Rp - Rs) / 2.0f, 1.0f - (Rs + Rp) / 2.0f, 0.0f),
			glm::vec3(0.0f, 0.0f, sqrt((1.0f - Rs) * (1.0f - Rp)))
		);
		res = Refl * res + Trans * polarInfoList[i].S;
	}
	glm::vec3 s_cam = glm::normalize(glm::cross(glm::vec3(0.0f, 0.0f, 1.0f), initDir));
	//phi = M_PI * 0.5f; // acos(0)
	phi = acos(glm::dot(s_cam, s)) * glm::sign(glm::dot(glm::cross(s, s_cam), initDir));

	return GetRotationMatrix(phi) * res;
}

__global__ void InitCuRand(int seed)
{
	const int x = threadIdx.x + blockIdx.x * blockDim.x;
	const int y = threadIdx.y + blockIdx.y * blockDim.y;
	if (x >= resX || y >= resY)
		return;

	curand_init(seed, x + y * resX, 0, &state[x + y * resX]);
}

void InitCUDA()
{
	float h_pi = glm::pi<float>();
	gpuErrchk(cudaMemcpyToSymbol(PI, &h_pi, sizeof(float)));
	gpuErrchk(cudaDeviceSetLimit(cudaLimitStackSize, 1024 * 8));
}

void CUDASetResolution(int x, int y)
{
	gpuErrchk(cudaMemcpyToSymbol(resX, &x, sizeof(unsigned)));
	gpuErrchk(cudaMemcpyToSymbol(resY, &y, sizeof(unsigned)));

	curandState_t* d_randState;
	gpuErrchk(cudaMalloc((void**)&d_randState, x * y * sizeof(curandState_t)));
	gpuErrchk(cudaMemcpyToSymbol(state, &d_randState, sizeof(d_randState)));

	srand(time(0));
	int seed = rand();
	dim3 blockDim(16, 16, 1), gridDim(x / blockDim.x + 1, y / blockDim.y + 1, 1);
	InitCuRand << < gridDim, blockDim >> > (seed);
	gpuErrchk(cudaGetLastError());
	gpuErrchk(cudaDeviceSynchronize());
}

void CUDASetTraceDepth(int depth)
{
	gpuErrchk(cudaMemcpyToSymbol(maxDepth, &depth, sizeof(unsigned)));
}

void CUDASetCamera(float* pos, float* dir, float* up)
{
	gpuErrchk(cudaMemcpyToSymbol(_camPos, pos, sizeof(float) * 3));
	gpuErrchk(cudaMemcpyToSymbol(_camDir, dir, sizeof(float) * 3));
	gpuErrchk(cudaMemcpyToSymbol(_camUp, up, sizeof(float) * 3));
}

void CUDASetProjection(float f, float fovy)
{
	gpuErrchk(cudaMemcpyToSymbol(camFocal, &f, sizeof(float)));
	gpuErrchk(cudaMemcpyToSymbol(camFovy, &fovy, sizeof(float)));
}

void CUDASetBVH(GPUBVHNode* nodes, int size)
{
	gpuErrchk(cudaMemcpyToSymbol(bvhSize, &size, sizeof(unsigned)));

	BVHNode* d_Nodes;
	gpuErrchk(cudaMalloc((void**)&d_Nodes, size * sizeof(GPUBVHNode)));
	gpuErrchk(cudaMemcpy(d_Nodes, nodes, size * sizeof(GPUBVHNode), cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpyToSymbol(bvh, &d_Nodes, sizeof(GPUBVHNode*)));
}

void CUDALoadTextures(const std::vector<Image*>& texVec)
{
	int size = texVec.size();
	numTextures = size;

	GPUImage* h_Imgs = new GPUImage[size];
	for (int i = 0; i < size; i++)
	{
		int w = texVec[i]->width();
		int h = texVec[i]->height();
		h_Imgs[i].width = w;
		h_Imgs[i].height = h;

		int size = w * h * 4 * sizeof(float);
		float* h_data = new float[size];
		memcpy(h_data, texVec[i]->data(), size);
		float* d_data;
		gpuErrchk(cudaMalloc(&d_data, size));
		gpuErrchk(cudaMemcpy(d_data, h_data, size, cudaMemcpyHostToDevice));
		delete[] h_data;
		h_Imgs->data = d_data;
	}

	GPUImage* d_Imgs;
	gpuErrchk(cudaMalloc(&d_Imgs, size * sizeof(GPUImage)));
	gpuErrchk(cudaMemcpy(d_Imgs, h_Imgs, size * sizeof(GPUImage), cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpyToSymbol(textures, &d_Imgs, size * sizeof(GPUImage*)));
	delete[] h_Imgs;
}

__device__ glm::vec4 CUDATex2D(const GPUImage& image, const glm::vec2& uv)
{
	if (uv.x > 1.0f || uv.x < 0.0f || uv.y > 1.0f || uv.y < 0.0f)
		return glm::vec4(0.0f);

	int w = image.width;
	int h = image.height;

	glm::ivec2 coord = glm::ivec2(w * uv.x, h * uv.y);
	float* p = image.data + (4 * (coord.y * w + coord.x));

	return glm::vec4(p[0], p[1], p[2], p[3]);
}

__device__ bool IsSameSide(const glm::vec3& p1, const glm::vec3& p2, const glm::vec3& a, const glm::vec3& b)
{
	glm::vec3 ba = b - a;
	glm::vec3 cp1 = glm::cross(ba, (p1 - a));
	glm::vec3 cp2 = glm::cross(ba, (p2 - a));

	return (glm::dot(cp1, cp2) >= 0);
}

__device__ bool IsInside(const glm::vec3& p, const glm::vec3& a, const glm::vec3& b, const glm::vec3& c)
{
	return (IsSameSide(p, a, b, c) && IsSameSide(p, b, a, c) && IsSameSide(p, c, a, b));
}

__device__ bool IntersectBox(const glm::vec3& ro, const glm::vec3& rd, const glm::vec3& bMin, const glm::vec3& bMax)
{
	glm::vec3 tMin = (bMin - ro) / rd;
	glm::vec3 tMax = (bMax - ro) / rd;
	glm::vec3 t1 = glm::min(tMin, tMax);
	glm::vec3 t2 = glm::max(tMin, tMax);
	float tNear = glm::max(glm::max(t1.x, t1.y), t1.z);
	float tFar = glm::min(glm::min(t2.x, t2.y), t2.z);
	if (tNear >= tFar)
		return false;
	return true;
}

__device__ bool Hit(const glm::vec3& ro, const glm::vec3& rd, Triangle& triangleOut, float& distOut)
{
	if (bvhSize <= 1)
		return false;

	bool res = false;

	distOut = float(0xFFFF);

	GPUBVHNode* stack[64];
	GPUBVHNode** pStack = stack;
	*pStack++ = NULL;

	GPUBVHNode* currNode = bvh;
	int stackIndex = 1;
	do
	{
		if (IntersectBox(ro, rd, currNode->box.min, currNode->box.max))
		{
			if (currNode->rightOffset == -1) // leaf
			{
				if (glm::dot(rd, currNode->triangle.normal) != 0.0f)
				{
					float d = glm::dot((currNode->triangle.v1 - ro), currNode->triangle.normal) / glm::dot(rd, currNode->triangle.normal);
					if (d >= 0)
					{
						glm::vec3 p = ro + rd * d;
						if (IsInside(p, currNode->triangle.v1, currNode->triangle.v2, currNode->triangle.v3))
						{
							if (d < distOut)
							{
								distOut = d;
								triangleOut = currNode->triangle;
							}
							res = true;
						}
					}
				}
				currNode = *--pStack;
				stackIndex--;
			}
			else // interier
			{
				GPUBVHNode* left = &(bvh[currNode->nodeIndex + 1]);
				GPUBVHNode* right = &(bvh[currNode->nodeIndex + currNode->rightOffset]);
				currNode = left;
				*pStack++ = right;
				stackIndex++;
			}
		}
		else
		{
			currNode = *--pStack;
			stackIndex--;
		}
	} while (stackIndex > 0 && stackIndex < 64);

	return res;
}

__device__ glm::vec2 GetUV(const glm::vec3& p, const Triangle& t)
{
	glm::vec3 v2 = p - t.v1;
	float d20 = glm::dot(v2, t.barycentricInfo.v0);
	float d21 = glm::dot(v2, t.barycentricInfo.v1);

	float alpha = (t.barycentricInfo.d11 * d20 - t.barycentricInfo.d01 * d21) *
		t.barycentricInfo.invDenom;
	float beta = (t.barycentricInfo.d00 * d21 - t.barycentricInfo.d01 * d20) *
		t.barycentricInfo.invDenom;

	return (1.0f - alpha - beta) * t.uv1 + alpha * t.uv2 + beta * t.uv3;
}

__device__ glm::vec3 GetSmoothNormal(const glm::vec3& p, const Triangle& t)
{
	glm::vec3 v2 = p - t.v1;
	float d20 = glm::dot(v2, t.barycentricInfo.v0);
	float d21 = glm::dot(v2, t.barycentricInfo.v1);

	float alpha = (t.barycentricInfo.d11 * d20 - t.barycentricInfo.d01 * d21) *
		t.barycentricInfo.invDenom;
	float beta = (t.barycentricInfo.d00 * d21 - t.barycentricInfo.d01 * d20) *
		t.barycentricInfo.invDenom;

	glm::vec3 n = (1.0f - alpha - beta) * t.n1 + alpha * t.n2 + beta * t.n3;
	glm::vec3 res = glm::normalize(glm::vec3(n.x, -n.y, n.z));
	return glm::normalize(n);
}

__device__ glm::vec3 reflect(glm::vec3 I, glm::vec3 N)
{
	return I - N * glm::dot(N, I) * glm::vec3(2);
}

__device__ glm::vec3 Trace(glm::vec3 ro, glm::vec3 rd, int& depth, bool& inside, curandState_t& state,
	GPUPolarInfo* polarInfoList)
{
	float d = 0.0f;
	Triangle t;
	if (Hit(ro, rd, t, d))
	{
		Material& mat = t.material;
		glm::vec3 p = ro + rd * d;
		glm::vec2 uv = GetUV(p, t);
		glm::vec3 n = t.normal;
		if (t.smoothing)
			n = GetSmoothNormal(p, t);
		if (glm::dot(n, rd) > 0.0f)
			n = -n;
		if (mat.normalTexId != -1)
		{
			glm::mat3 TBN = glm::mat3(t.tangent, t.bitangent, n);
			glm::vec3 nt = glm::vec3(CUDATex2D(textures[mat.normalTexId], uv)) * 2.0f - 1.0f;
			if (nt.z < 0.0f)
				nt = glm::vec3(nt.x, nt.y, 0.0f);
			nt = glm::normalize(nt);
			n = glm::normalize(TBN * nt);
		}
		p += n * EPSILON;

		if (depth < maxDepth * 2)
		{
			depth++;
			// Russian Roulette Path Termination
			float prob = glm::min(0.95f, glm::max(glm::max(mat.baseColor.x, mat.baseColor.y), mat.baseColor.z));
			if (depth >= maxDepth)
			{
				if (fabs(curand_uniform(&state)) > prob)
					return mat.emissive * mat.emissiveIntensity;
			}

			glm::vec3 r = reflect(rd, n);
			glm::vec3 reflectDir = r;
			if (mat.type == MaterialType::SPECULAR)
				reflectDir = r;
			else if (mat.type == MaterialType::DIFFUSE)
			{
				// Monte Carlo Integration
				glm::vec3 u = glm::abs(n.x) < 1.0f - EPSILON ? glm::cross(glm::vec3(1.0f, 0.0f, 0.0f), n) : glm::cross(glm::vec3(1.0f), n);
				u = glm::normalize(u);
				glm::vec3 v = glm::normalize(glm::cross(u, n));
				float w = curand_uniform(&state), theta = curand_uniform(&state);
				// uniformly sampling on hemisphere
				reflectDir = w * cosf(2.0f * PI * theta) * u + w * sinf(2.0f * PI * theta) * v + glm::sqrt(1.0f - w * w) * n;
				reflectDir = glm::normalize(reflectDir);
			}
			else if (mat.type == MaterialType::GLOSSY)
			{
				// Monte Carlo Integration
				glm::vec3 u = fabs(n.x) < 1 - FLT_EPSILON ? glm::cross(glm::vec3(1, 0, 0), r) : glm::cross(glm::vec3(1), r);
				u = glm::normalize(u);
				glm::vec3 v = glm::cross(u, r);
				float w = curand_uniform(&state) * mat.roughness, theta = curand_uniform(&state);
				// wighted sampling on hemisphere
				reflectDir = w * cosf(2 * PI * theta) * u + w * sinf(2 * PI * theta) * v + sqrtf(1 - w * w) * r;
			}
			else if (mat.type == MaterialType::GLASS)
			{
				float nc = 1.0f, ng = 1.5f;
				// Snells law
				float eta = inside ? ng / nc : nc / ng;
				float r0 = (nc - ng) / (nc + ng);
				r0 = r0 * r0;
				float c = fabs(glm::dot(rd, n));
				float k = 1.0f - eta * eta * (1.0f - c * c);
				if (k < 0.0f)
					reflectDir = r;
				else
				{
					// Shilick's approximation of Fresnel's equation
					float re = r0 + (1.0f - r0) * (1.0f - c) * (1.0f - c);
					if (fabs(curand_uniform(&state)) < re)
						reflectDir = r;
					else
					{
						reflectDir = glm::normalize(eta * rd - (eta * glm::dot(n, rd) + sqrtf(k)) * n);
						p -= n * EPSILON * 2.0f;
						inside = !inside;
					}
				}
			}

			polarInfoList[depth - 1].indir = rd;
			polarInfoList[depth - 1].outdir = reflectDir;
			polarInfoList[depth - 1].normal = glm::normalize(reflectDir - rd);
			if (fabs(1.0f - glm::dot(reflectDir, rd)) < EPSILON)
				polarInfoList[depth - 1].normal = n;
			polarInfoList[depth - 1].R = GetRsRp(1.0f, mat.ior, glm::dot(polarInfoList[depth - 1].normal, rd));
			float theta = 2.0f * PI * curand_uniform(&state);
			float intensity = mat.intensity;
			if (mat.intensityTexId != -1)
			{
				// TODO map intensity value here
				intensity = CUDATex2D(textures[mat.intensityTexId], uv).x;
			}
			polarInfoList[depth - 1].S = glm::vec3(1.0f, cosf(theta), sinf(theta)) * intensity;
			
			return mat.emissive * mat.emissiveIntensity + Trace(p, reflectDir, depth, inside, state, polarInfoList) * mat.baseColor;
		}
	}

	return glm::vec3(0.0f);
}

__global__ void RenderPixel(float* img, float* data)
{
	const int x = threadIdx.x + blockIdx.x * blockDim.x;
	const int y = threadIdx.y + blockIdx.y * blockDim.y;
	if (x >= resX || y >= resY)
		return;
	const int index = x + (resY - y - 1)*resX;
	curandState_t localState = state[index];

	// Position world space image plane
	glm::vec3 camPos = glm::vec3(_camPos[0], _camPos[1], _camPos[2]);
	glm::vec3 camDir = glm::vec3(_camDir[0], _camDir[1], _camDir[2]);
	glm::vec3 camUp = glm::vec3(_camUp[0], _camUp[1], _camUp[2]);
	glm::vec3 imgCenter = camPos + camDir * camFocal;
	float imgHeight = 2.0f * camFocal * tan((camFovy / 2.0f) * PI / 180.0f);
	float aspect = (float)resX / (float)resY;
	float imgWidth = imgHeight * aspect;
	float deltaX = imgWidth / (float)resX;
	float deltaY = imgHeight / (float)resY;
	glm::vec3 camRight = glm::normalize(glm::cross(camUp, camDir));

	// Starting at top left
	glm::vec3 topLeft = imgCenter - camRight * (imgWidth * 0.5f);
	topLeft += camUp * (imgHeight * 0.5f);
	glm::vec3 pixel = topLeft - camUp * (float(y) * deltaY) + camRight * (float(x) * deltaX);

	/* ----- GATHER POLAR INFO ----- */
	GPUPolarInfo polarInfoList[10];
	/* ----- GATHER POLAR INFO ----- */

	glm::vec3 rayDir = glm::normalize(pixel - camPos);
	int depth = 0;
	bool inside = false;
	glm::vec3 color = Trace(camPos, rayDir, depth, inside, localState, polarInfoList);

	/* ----- CALCULAT POLAR RESULT ----- */
	glm::vec3 polarResult = CalculatePolarResult(rayDir, polarInfoList, depth, localState);
	if (isnan(polarResult.x) || isnan(polarResult.y) || isnan(polarResult.z))
		polarResult = glm::vec3(0.0f);
	color = polarResult;
	/* ----- CALCULAT POLAR RESULT ----- */

	// Draw
	glm::vec3 colorR = color * 5.0f;
	glm::vec3 colorGB = 1.0f - (color + 0.025f) * 20.0f;
	color = glm::vec3(colorR.x, colorGB.y, colorGB.z);

	glm::vec3 preColor = glm::vec3();
	memcpy(&preColor[0], img + 3 * index, 3 * sizeof(float));
	color = (preColor * float(samples - 1) + color) / float(samples);
	color = glm::clamp(color, glm::vec3(0.0f), glm::vec3(1.0f));
	memcpy(img + 3 * index, &color[0], 3 * sizeof(float));

	// Data
	glm::vec3 preRes = glm::vec3();
	memcpy(&preRes[0], data + 3 * index, 3 * sizeof(float));
	polarResult = (preRes * float(samples - 1) + polarResult) / float(samples);
	memcpy(data + 3 * index, &polarResult[0], 3 * sizeof(float));

	state[index] = localState;
}

void CUDARenderFrame(int w, int h, float* img, float* data, int& h_samples)
{
	gpuErrchk(cudaMemcpyToSymbol(samples, &h_samples, sizeof(int)));
	dim3 blockDim(16, 16, 1), gridDim(w / blockDim.x + 1, h / blockDim.y + 1, 1);
	RenderPixel << <gridDim, blockDim >> > (img, data);
	gpuErrchk(cudaGetLastError());
	gpuErrchk(cudaDeviceSynchronize());
}

void CUDAReset()
{
	if (bvh != 0)
	{
		gpuErrchk(cudaFree(bvh));
		bvh = 0;
	}

	if (textures != 0)
	{
		for (int i = 0; i < numTextures; i++)
			gpuErrchk(cudaFree(textures[i].data));
		gpuErrchk(cudaFree(textures));
		textures = 0;
		numTextures = 0;
	}
}

void CUDAFinish()
{
	gpuErrchk(cudaFree(state));
}
