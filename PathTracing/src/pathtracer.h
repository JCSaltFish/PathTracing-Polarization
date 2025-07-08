#ifndef __PATHTRACER_H__
#define __PATHTRACER_H__

#include <string>
#include <vector>
#include <random>

#include <GL/glew.h>
#include <glm/glm.hpp>

#include "mesh.h"

enum class MaterialType
{
	DIFFUSE,
	SPECULAR,
	GLOSSY,
	GLASS
};

class IntensityData
{
public:
	IntensityData(const std::string& file);
	~IntensityData() {}

	float Read(const glm::vec2& uv)
	{
		if (uv.x > 1.0f || uv.x < 0.0f || uv.y > 1.0f || uv.y < 0.0f)
			return 0.0f;
		glm::ivec2 coord = glm::ivec2(mWidth * uv.x, mHeight * uv.y);
		return mData[coord.y * mWidth + coord.x];
	}

private:
	std::vector<float> mData;
	int mWidth;
	int mHeight;
};

struct Material
{
	MaterialType type = MaterialType::DIFFUSE;
	glm::vec3 baseColor = glm::vec3(1.0f);
	float roughness = 0.0f;
	glm::vec3 emissive = glm::vec3(0.0f);

	int normalTexId = -1;

	float ior = 1.0f;

	float intensity = 0.0f;
	int intensityTexId = -1;
	//粗糙度后加
	int roughnessTexId = -1;

	IntensityData* pIntensityData = nullptr;

	~Material()
	{
		if (pIntensityData)
			delete pIntensityData;
	}
};

namespace PathTracerLoader
{
	struct Element
	{
		std::string name;
		Material material;

		Element()
		{
			name = "";
		}

		Element(const std::string& name)
		{
			this->name = name;
		}
	};

	struct Object
	{
		std::string name;
		std::vector<Element> elements;

		Object()
		{
			name = "";
		}

		Object(const std::string& name)
		{
			this->name = name;
		}
	};
}

struct PolarInfo
{
	glm::vec3 normal;
	glm::vec3 indir;
	glm::vec3 outdir;

	glm::vec2 R;
	glm::vec3 S;
};

class PathTracer
{
private:
	std::vector<Triangle> mTriangles;
	BVHNode* mBvh;

	std::vector<PathTracerLoader::Object> mLoadedObjects;
	std::vector<Image*> mLoadedTextures;

	glm::ivec2 mResolution;
	GLubyte* mOutImg;
	float* mTotalImg;
	int mMaxDepth;

	glm::vec3 mCamPos;
	glm::vec3 mCamDir;
	glm::vec3 mCamUp;
	float mCamFocal;
	float mCamFovy;

	int mSamples;
	bool mNeedReset;
	bool mExit;

	std::mt19937 mRng;

	float* mPolarData;

public:
	PathTracer();
	~PathTracer();

private:
	/* ----- POLAR FUNCTIONS ----- */
	const glm::vec2 GetRsRp(float n1, float n2, float cos_i) const;
	const glm::mat3 GetRotationMatrix(float phi) const;
	const glm::vec3 CalculatePolarResult(const glm::vec3& initDir, const std::vector<PolarInfo*>& polarInfoList);
	/* ----- POLAR FUNCTIONS ----- */

	const glm::vec3 IntersectTriangle
	(
		const glm::vec3& ro, const glm::vec3& rd,
		const glm::vec3& v0, const glm::vec3& v1, const glm::vec3& v2
	) const;
	const bool LinearHit(const glm::vec3& ro, const glm::vec3& rd, Triangle& triangleOut, float& distOut);

	const float Rand();
	const glm::vec2 GetUV(const glm::vec3& p, const Triangle& t) const;
	const glm::vec3 GetSmoothNormal(const glm::vec3& p, const Triangle& t) const;
	const glm::vec3 Trace(const glm::vec3& ro, const glm::vec3& rd, std::vector<PolarInfo*>& polarInfoList, int depth = 0, bool inside = false);

public:
	void LoadObject(const std::string& file, const glm::mat4& model);
	void SetNormalTextureForElement(int objId, int elementId, const std::string& file);
	//粗糙度后加
	void SetRoughnessTextureForElement(int objId, int elementId, const std::string& file);
	//
	void SetMaterial(int objId, int elementId, Material& material);
	void BuildBVH();
	void ResetImage();
	void ClearScene();

	const int GetSamples() const;
	const int GetTriangleCount() const;
	const int GetTraceDepth() const;
	void SetTraceDepth(int depth);
	void SetOutImage(GLubyte* out);
	void SetResolution(const glm::ivec2& res);
	const glm::ivec2 GetResolution() const;
	std::vector<PathTracerLoader::Object> GetLoadedObjects() const;

	void SetPolarData(float* data);
	void SetIntensityTextureForElement(int objId, int elementId, const std::string& file);
	void SetIntensityDataForElement(int objId, int elementId, const std::string& file);

	void SetCamera(const glm::vec3& pos, const glm::vec3& dir, const glm::vec3& up);
	void SetProjection(float f, float fovy);
	void RenderFrame();
	void Exit();
};

#endif
