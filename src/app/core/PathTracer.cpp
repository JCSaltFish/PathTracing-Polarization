/**
 * @file PathTracer.cpp
 * @brief Implementation file for the PathTracer class.
 */

#include "app/core/PathTracer.h"

#include "app/AppTextureManager.h"
#include "utils/Logger.hpp"
#include "utils/Flags.hpp"
#include "res/ShaderStringsUtils.hpp"

int PathTracer::init() {
    if (!m_renderer) {
        Logger() << "Invalid renderer in PathTracer::init";
        return 1;
    }

    // Load shader
    try {
        m_computeShader = m_renderer->createShader
        (
            GfxShaderStage::COMPUTE,
            ShaderStrings::get("pathTracer.comp")
        );
    } catch (GfxShaderException& e) {
        Logger() << "Shader compilation error in PathTracer::init: " << e.what();
        return 1;
    }

    /* Initialize descriptors and UBOs */
    m_descriptors.u_outImage.binding = 0;
    m_descriptors.u_outImage.type = GfxDescriptorType::STORAGE_IMAGE;
    m_descriptors.u_outImage.stages.set(GfxShaderStage::COMPUTE);

    m_descriptors.u_scene.binding = 1;
    m_descriptors.u_scene.type = GfxDescriptorType::UNIFORM_BUFFER;
    m_descriptors.u_scene.stages.set(GfxShaderStage::COMPUTE);
    m_uboScene = m_renderer->createBuffer(
        sizeof(UScene),
        GfxBufferUsage::UNIFORM_BUFFER,
        GfxBufferProp::DYNAMIC
    );
    if (!m_uboScene) {
        Logger() << "Failed to create scene UBO in PathTracer::init";
        return 1;
    }

    m_descriptors.u_camera.binding = 2;
    m_descriptors.u_camera.type = GfxDescriptorType::UNIFORM_BUFFER;
    m_descriptors.u_camera.stages.set(GfxShaderStage::COMPUTE);
    m_uboCamera = m_renderer->createBuffer(
        sizeof(UCamera),
        GfxBufferUsage::UNIFORM_BUFFER,
        GfxBufferProp::DYNAMIC
    );
    if (!m_uboCamera) {
        Logger() << "Failed to create camera UBO in PathTracer::init";
        return 1;
    }

    m_descriptors.u_textures.binding = 3;
    m_descriptors.u_textures.type = GfxDescriptorType::SAMPLERS;
    m_descriptors.u_textures.stages.set(GfxShaderStage::COMPUTE);

    m_descriptors.b_vertices.binding = 4;
    m_descriptors.b_vertices.type = GfxDescriptorType::STORAGE_BUFFER;
    m_descriptors.b_vertices.stages.set(GfxShaderStage::COMPUTE);

    m_descriptors.b_triangles.binding = 5;
    m_descriptors.b_triangles.type = GfxDescriptorType::STORAGE_BUFFER;
    m_descriptors.b_triangles.stages.set(GfxShaderStage::COMPUTE);

    m_descriptors.b_materials.binding = 6;
    m_descriptors.b_materials.type = GfxDescriptorType::STORAGE_BUFFER;
    m_descriptors.b_materials.stages.set(GfxShaderStage::COMPUTE);

    m_descriptors.b_BVH.binding = 7;
    m_descriptors.b_BVH.type = GfxDescriptorType::STORAGE_BUFFER;
    m_descriptors.b_BVH.stages.set(GfxShaderStage::COMPUTE);

    return 0;
}

void PathTracer::term() {
    if (!m_renderer)
        return;

    clearScene();

    m_renderer->destroyBuffer(m_uboScene);
    m_renderer->destroyBuffer(m_uboCamera);
    m_renderer->destroyShader(m_computeShader);

    m_descriptors = {};
}

int PathTracer::buildScene(const DbObjHandle& hScene) {
    if (!hScene.isValid() || hScene.getType() != PtScene::TYPE_NAME) {
        Logger() << "Invalid scene handle in PathTracer::buildScene";
        return 1;
    }
    if (!m_renderer) {
        Logger() << "Invalid renderer in PathTracer::buildScene";
        return 1;
    }

    PtScene::getResolution(hScene, m_resolutionX, m_resolutionY);
    BufferData bufferData = {};

    /* Load models */
    loadModels(hScene, bufferData);

    m_renderer->waitDeviceIdle();

    /* Create pipeline */
    m_descriptors.u_textures.size = static_cast<int>(bufferData.textures.size());
    if (m_pipeline)
        m_renderer->destroyPipeline(m_pipeline);
    m_pipeline = m_renderer->createPipeline(
        { m_computeShader },
        {
            {
                m_descriptors.u_outImage,
                m_descriptors.u_scene,
                m_descriptors.u_camera,
                m_descriptors.u_textures,
                m_descriptors.b_vertices,
                m_descriptors.b_triangles,
                m_descriptors.b_materials,
                m_descriptors.b_BVH
            }
        }
    );
    if (!m_pipeline) {
        Logger() << "Failed to create pipeline in PathTracer::buildScene";
        return 1;
    }

    /* Create GPU buffers */
    if (createBuffers(bufferData)) {
        Logger() << "Failed to create GPU buffers in PathTracer::buildScene";
        return 1;
    }

    /* Create output and display image */
    if (m_outImage)
        m_renderer->destroyImage(m_outImage);
    GfxImageInfo outImgInfo = {};
    outImgInfo.width = m_resolutionX;
    outImgInfo.height = m_resolutionY;
    outImgInfo.format = GfxFormat::R32G32B32A32_SFLOAT;
    outImgInfo.usages = GfxImageUsage::STORAGE_IMAGE;
    m_outImage = m_renderer->createImage(outImgInfo);
    if (!m_outImage) {
        Logger() << "Failed to create output image in PathTracer::buildScene";
        return 1;
    }
    if (m_dspImageFront)
        m_renderer->destroyImage(m_dspImageFront);
    if (m_dspImageBack)
        m_renderer->destroyImage(m_dspImageBack);
    GfxImageInfo dspImgInfo = {};
    dspImgInfo.width = m_resolutionX;
    dspImgInfo.height = m_resolutionY;
    dspImgInfo.format = GfxFormat::R32G32B32A32_SFLOAT;
    dspImgInfo.usages.set(GfxImageUsage::COLOR_ATTACHMENT);
    dspImgInfo.usages.set(GfxImageUsage::SAMPLED_TEXTURE);
    m_dspImageFront = m_renderer->createImage(dspImgInfo);
    m_dspImageBack = m_renderer->createImage(dspImgInfo);
    if (!m_dspImageFront || !m_dspImageBack) {
        Logger() << "Failed to create display image in PathTracer::buildScene";
        return 1;
    }

    /* Create descriptor set binding */
    if (m_descriptorSetBinding)
        m_renderer->destroyDescriptorSetBinding(m_descriptorSetBinding);
    std::vector<GfxDescriptorBinding> bindings = {};
    bindings.reserve(8);
    bindings.push_back({ m_descriptors.u_outImage, m_outImage });
    bindings.push_back({ m_descriptors.u_scene, m_uboScene });
    bindings.push_back({ m_descriptors.u_camera, m_uboCamera });
    bindings.push_back({ m_descriptors.u_textures, bufferData.textures });
    bindings.push_back({ m_descriptors.b_vertices, m_ssboVertex });
    bindings.push_back({ m_descriptors.b_triangles, m_ssboTriangle });
    bindings.push_back({ m_descriptors.b_materials, m_ssboMaterial });
    bindings.push_back({ m_descriptors.b_BVH, m_ssboBVH });
    m_descriptorSetBinding = m_renderer->createDescriptorSetBinding(m_pipeline, 0, bindings);

    /* Load scene settings and update UBOs */
    UScene u_scene = {};
    u_scene.resX = m_resolutionX;
    u_scene.resY = m_resolutionY;
    u_scene.traceDepth = PtScene::getTraceDepth(hScene);
    m_currentSample = 0;
    if (m_renderer->updateBufferData(m_uboScene, 0, sizeof(u_scene), &u_scene)) {
        Logger() << "Failed to update scene UBO in PathTracer::buildScene";
        return 1;
    }

    using namespace Math;
    UCamera u_camera = {};
    PtScene::Camera sceneCam = PtScene::getCamera(hScene);
    u_camera.pos = Vec4(sceneCam.position, 1.0f);
    Mat4 rotX = rotate(Mat4(1.0f), sceneCam.rotation.x, Vec3(1.0f, 0.0f, 0.0f));
    Mat4 rotY = rotate(Mat4(1.0f), sceneCam.rotation.y, Vec3(0.0f, 1.0f, 0.0f));
    Mat4 rotZ = rotate(Mat4(1.0f), sceneCam.rotation.z, Vec3(0.0f, 0.0f, 1.0f));
    Mat4 rot = rotZ * rotY * rotX;
    u_camera.dir = rot * Vec4(0.0f, 0.0f, 1.0f, 0.0f);
    u_camera.up = rot * Vec4(0.0f, 1.0f, 0.0f, 0.0f);
    u_camera.focusDist = sceneCam.focusDist;
    u_camera.fStop = sceneCam.fStop;
    if (m_renderer->updateBufferData(m_uboCamera, 0, sizeof(u_camera), &u_camera)) {
        Logger() << "Failed to update camera UBO in PathTracer::buildScene";
        return 1;
    }

    return 0;
}

void PathTracer::clearScene() {
    if (!m_renderer)
        return;
    m_renderer->waitDeviceIdle();

    if (m_descriptorSetBinding) {
        m_renderer->destroyDescriptorSetBinding(m_descriptorSetBinding);
        m_descriptorSetBinding = nullptr;
    }
    if (m_outImage) {
        m_renderer->destroyImage(m_outImage);
        m_outImage = nullptr;
    }
    if (m_dspImageFront) {
        m_renderer->destroyImage(m_dspImageFront);
        m_dspImageFront = nullptr;
    }
    if (m_dspImageBack) {
        m_renderer->destroyImage(m_dspImageBack);
        m_dspImageBack = nullptr;
    }

    if (m_ssboVertex) {
        m_renderer->destroyBuffer(m_ssboVertex);
        m_ssboVertex = nullptr;
    }
    if (m_ssboTriangle) {
        m_renderer->destroyBuffer(m_ssboTriangle);
        m_ssboTriangle = nullptr;
    }
    if (m_ssboMaterial) {
        m_renderer->destroyBuffer(m_ssboMaterial);
        m_ssboMaterial = nullptr;
    }
    if (m_ssboBVH) {
        m_renderer->destroyBuffer(m_ssboBVH);
        m_ssboBVH = nullptr;
    }

    if (m_pipeline) {
        m_renderer->destroyPipeline(m_pipeline);
        m_pipeline = nullptr;
    }

    m_currentSample = 0;
}

int PathTracer::renderFrame() {
    m_renderer->bindPipeline(m_pipeline);

    // Update current sample in UBO
    m_currentSample++;
    int err = m_renderer->updateBufferData(
        m_uboScene,
        static_cast<int>(offsetof(UScene, currentSample)),
        static_cast<int>(sizeof(uint32_t)),
        &m_currentSample
    );
    if (err)
        return 1;
    m_renderer->bindDescriptorSetBinding(m_descriptorSetBinding);

    // Dispatch compute shader
    m_renderer->dispatchCompute(
        static_cast<int>(std::ceil(static_cast<float>(m_resolutionX) / 32.0f)),
        static_cast<int>(std::ceil(static_cast<float>(m_resolutionY) / 32.0f)),
        1
    );
    m_renderer->memoryBarrier();

    // Copy output image to display image
    m_renderer->copyImage(m_outImage, m_dspImageFront, m_resolutionX, m_resolutionY);

    return 0;
}

uint32_t PathTracer::getCurrentSample() const {
    return m_currentSample;
}

GfxImage PathTracer::getCurrentDisplayImage() const {
    if (m_dspImageFront)
        return m_dspImageFront;
    else
        return AppTextureManager::instance().getDefaultTexture();
}

std::array<GfxImage, 2> PathTracer::getDisplayImages() const {
    return { m_dspImageFront, m_dspImageBack };
}

void PathTracer::markDisplayImageReady() {
    m_dspImgSwapPending.store(true, std::memory_order_release);
}

void PathTracer::syncDisplayImage() {
    if (m_dspImgSwapPending.load(std::memory_order_acquire)) {
        std::swap(m_dspImageFront, m_dspImageBack);
        m_dspImgSwapPending.store(false, std::memory_order_release);
    }
}

int PathTracer::getImageData(
    std::vector<unsigned char>& pixels,
    int& width,
    int& height
) const {
    if (!m_renderer || !m_outImage)
        return 1;
    size_t pixelCount =
        static_cast<size_t>(m_resolutionX) * static_cast<size_t>(m_resolutionY);
    std::vector<float> data(pixelCount * 4);
    pixels.resize(pixelCount * 4 * sizeof(unsigned char));
    if (m_renderer->getImageData(m_outImage, data.data()))
        return 1;
    for (size_t i = 0; i < pixelCount * 4; i++) {
        float value = data[i];
        value = std::clamp(value, 0.0f, 1.0f);
        pixels[i] = static_cast<unsigned char>(value * 255.0f);
    }
    width = m_resolutionX;
    height = m_resolutionY;
    return 0;
}

void PathTracer::render() {
    m_rendering = true;
}

void PathTracer::pause() {
    m_rendering = false;
}

void PathTracer::stop() {
    m_rendering = false;
    m_currentSample = 0;
}

void PathTracer::restart() {
    m_currentSample = 0;
    m_rendering = true;
}

bool PathTracer::isRendering() const {
    return m_rendering;
}

bool PathTracer::isPaused() const {
    return !m_rendering && m_currentSample > 0;
}

void PathTracer::renderFinishCallback() {
    if (!m_renderFinishCb)
        return;
    m_renderFinishCb();
    m_renderFinishCb = nullptr;
}

void PathTracer::setRenderFinishCallback(const std::function<void()>& cb) {
    m_renderFinishCb = cb;
}

void PathTracer::loadModels(const DbObjHandle& hScene, BufferData& data) {
    for (const auto& hModel : PtScene::getModels(hScene)) {
        /* Load model data from file */
        std::string filename = PtModel::getFilePath(hModel);
        if (filename.empty()) {
            Logger() << "Model file path is empty for model ID: " << hModel.getID();
            continue;
        }
        Mesh::Model modelData = {};
        if (MeshLoader::loadOBJ(filename, modelData)) {
            Logger() << "Failed to load model file: " << filename;
            continue;
        }

        std::vector<DbObjHandle> meshHandles = PtModel::getMeshes(hModel);

        /* Pre-check mesh count */
        int meshCount = 0;
        for (const auto& meshData : modelData.meshes) {
            for (const auto& submeshData : meshData.submeshes)
                meshCount++;
        }
        if (meshHandles.size() != meshCount) {
            Logger() << "Mesh count mismatch for model: " << filename;
            continue;
        }

        /* Sync materials */
        uint32_t idxMaterial = static_cast<uint32_t>(data.materials.size());
        for (const auto& hMesh : meshHandles) {
            data.materials.push_back({});

            if (!hMesh.isValid() || hMesh.getType() != PtMesh::TYPE_NAME) {
                Logger() << "Invalid mesh handle in model: " << filename;
                continue;
            }
            DbObjHandle hMaterial = PtMesh::getMaterial(hMesh);
            if (!hMaterial.isValid() || hMaterial.getType() != PtMaterial::TYPE_NAME) {
                Logger() << "Invalid material handle in mesh ID: " << hMesh.getID();
                continue;
            }

            Material& material = data.materials.back();
            material.type = static_cast<int>(PtMaterial::getType(hMaterial));
            material.roughness = PtMaterial::getRoughness(hMaterial);
            material.intensity = PtMaterial::getIntensity(hMaterial);
            material.ior = PtMaterial::getIOR(hMaterial);
            material.flags = PtMaterial::getFlags(hMaterial).getValue();

            Flags<PtMaterial::MaterialFlag> materialFlags = material.flags;
            // Load textures
            auto getTexture = [](const std::string& texPath, std::vector<GfxImage>& textures) {
                GfxImage texture = AppTextureManager::instance().getTexture(texPath);
                if (texture != nullptr) {
                    textures.push_back(texture);
                    return static_cast<uint32_t>(textures.size() - 1);
                } else
                    return static_cast<uint32_t>(0);
                };
            if (materialFlags.check(PtMaterial::MaterialFlag::NORMAL_MAP)) {
                material.idxNormalTex =
                    getTexture(PtMaterial::getNormalTexPath(hMaterial), data.textures);
            }
            if (materialFlags.check(PtMaterial::MaterialFlag::ROUGHNESS_MAP)) {
                material.idxRoughnessTex =
                    getTexture(PtMaterial::getRoughnessTexPath(hMaterial), data.textures);
            }
            if (materialFlags.check(PtMaterial::MaterialFlag::INTENSITY_MAP)) {
                std::string texPath = PtMaterial::getIntensityTexPath(hMaterial);
                GfxImage texture = AppTextureManager::instance().getIntensityTexture(texPath);
                if (texture != nullptr) {
                    data.textures.push_back(texture);
                    material.idxIntensityTex = data.textures.size() - 1;
                }
            }
        }

        /* Pre-process model transformation */
        using namespace Math;
        Vec3 location = PtModel::getLocation(hModel);
        location = Vec3(-location.x, location.y, location.z);
        Mat4 t = translate(Mat4(1.0f), location);
        Vec3 rotation = PtModel::getRotation(hModel);
        rotation = Vec3(rotation.x, -rotation.y, -rotation.z);
        Mat4 rx = rotate(Mat4(1.0f), rotation.x, Vec3(1.0f, 0.0f, 0.0f));
        Mat4 ry = rotate(Mat4(1.0f), rotation.y, Vec3(0.0f, 1.0f, 0.0f));
        Mat4 rz = rotate(Mat4(1.0f), rotation.z, Vec3(0.0f, 0.0f, 1.0f));
        Mat4 s = scale(Mat4(1.0f), PtModel::getScale(hModel));
        Mat4 xform = t * rz * ry * rx * s;

        /* Process model data */
        for (int i = 0; i < modelData.meshes.size(); i++) {
            /* Process mesh data */
            const Mesh::Mesh& meshData = modelData.meshes[i];
            uint32_t vtxIdxOffset = static_cast<uint32_t>(data.vertices.size());

            // Vertices
            for (const auto& vtxData : meshData.vertices) {
                Vertex vtx = {};
                vtx.pos = xform * Vec4(vtxData.pos, 1.0f);
                vtx.normal = xform * Vec4(vtxData.normal, 0.0f);
                vtx.tangent = xform * Vec4(vtxData.tangent, 0.0f);
                vtx.texCoord = vtxData.texCoord;
                data.vertices.push_back(vtx);
            }

            // Triangles
            for (const auto& submeshData : meshData.submeshes) {
                if (submeshData.indices.size() < 3) {
                    idxMaterial++;
                    continue;
                }
                for (size_t i = 0; i < submeshData.indices.size() - 2; i += 3) {
                    Triangle t = {};
                    t.v0 = vtxIdxOffset + submeshData.indices[i + 0];
                    t.v1 = vtxIdxOffset + submeshData.indices[i + 1];
                    t.v2 = vtxIdxOffset + submeshData.indices[i + 2];
                    t.idxMaterial = idxMaterial;
                    data.triangles.push_back(t);
                }
                idxMaterial++;
            }
        }
    }

    /* Build scene BVH */
    BvhBuilder bvhBuilder;
    std::shared_ptr<BvhNode> bvh = bvhBuilder.build(data.vertices, data.triangles);
    BvhBufferizer bvhBufferizer;
    data.bvhBufferData = bvhBufferizer.bufferize(bvh.get());

    if (data.textures.empty())
        data.textures.push_back(AppTextureManager::instance().getDefaultTexture());
}

int PathTracer::createBuffers(const BufferData& data) {
    int err = 0;

    // Vertex buffer
    if (m_ssboVertex)
        m_renderer->destroyBuffer(m_ssboVertex);
    m_ssboVertex = m_renderer->createBuffer(
        static_cast<int>(sizeof(Vertex) * data.vertices.size()),
        GfxBufferUsage::STORAGE_BUFFER,
        GfxBufferProp::DYNAMIC
    );
    if (!m_ssboVertex)
        return 1;
    err = m_renderer->setBufferData(
        m_ssboVertex,
        static_cast<int>(sizeof(Vertex) * data.vertices.size()),
        data.vertices.data()
    );
    if (err)
        return 1;

    // Triangle buffer
    if (m_ssboTriangle)
        m_renderer->destroyBuffer(m_ssboTriangle);
    m_ssboTriangle = m_renderer->createBuffer(
        static_cast<int>(sizeof(Triangle) * data.triangles.size()),
        GfxBufferUsage::STORAGE_BUFFER,
        GfxBufferProp::DYNAMIC
    );
    if (!m_ssboTriangle)
        return 1;
    err = m_renderer->setBufferData(
        m_ssboTriangle,
        static_cast<int>(sizeof(Triangle) * data.triangles.size()),
        data.triangles.data()
    );
    if (err)
        return 1;

    // Material buffer
    if (m_ssboMaterial)
        m_renderer->destroyBuffer(m_ssboMaterial);
    m_ssboMaterial = m_renderer->createBuffer(
        static_cast<int>(sizeof(Material) * data.materials.size()),
        GfxBufferUsage::STORAGE_BUFFER,
        GfxBufferProp::DYNAMIC
    );
    if (!m_ssboMaterial)
        return 1;
    err = m_renderer->setBufferData(
        m_ssboMaterial,
        static_cast<int>(sizeof(Material) * data.materials.size()),
        data.materials.data()
    );
    if (err)
        return 1;

    // BVH buffer
    if (m_ssboBVH)
        m_renderer->destroyBuffer(m_ssboBVH);
    m_ssboBVH = m_renderer->createBuffer(
        static_cast<int>(sizeof(BufferBvhNode) * data.bvhBufferData.size()),
        GfxBufferUsage::STORAGE_BUFFER,
        GfxBufferProp::DYNAMIC
    );
    if (!m_ssboBVH)
        return 1;
    err = m_renderer->setBufferData(
        m_ssboBVH,
        static_cast<int>(sizeof(BufferBvhNode) * data.bvhBufferData.size()),
        data.bvhBufferData.data()
    );
    if (err)
        return 1;

    return 0;
}

std::unique_ptr<PathTracer::BvhNode> PathTracer::BvhBuilder::build
(
    const std::vector<Vertex>& vertices,
    const std::vector<Triangle>& triangles
) {
    m_triList.resize(triangles.size());
    m_triAABBs.resize(triangles.size());
    for (int i = 0; i < triangles.size(); i++) {
        m_triList[i] = i;
        m_triAABBs[i].merge(Math::Vec3(vertices[triangles[i].v0].pos));
        m_triAABBs[i].merge(Math::Vec3(vertices[triangles[i].v1].pos));
        m_triAABBs[i].merge(Math::Vec3(vertices[triangles[i].v2].pos));
        m_triAABBs[i].validate();
    }
    std::unique_ptr<BvhNode> root = std::make_unique<BvhNode>();
    buildRecursive(root.get(), 0, triangles.size());
    return root;
}

void PathTracer::BvhBuilder::buildRecursive(BvhNode* node, size_t triListOffset, size_t triCount) {
    for (int i = triListOffset; i < triListOffset + triCount; i++)
        node->aabb.merge(m_triAABBs[m_triList[i]]);

    /* Build leaves */
    if (triCount == 0)
        return;
    else if (triCount == 1) {
        node->left = std::make_unique<BvhNode>();
        node->left->aabb = m_triAABBs[m_triList[triListOffset]];
        node->left->idxTriangle = m_triList[triListOffset];
        return;
    } else if (triCount == 2) {
        node->left = std::make_unique<BvhNode>();
        node->left->aabb = m_triAABBs[m_triList[triListOffset + 0]];
        node->left->idxTriangle = m_triList[triListOffset + 0];
        node->right = std::make_unique<BvhNode>();
        node->right->aabb = m_triAABBs[m_triList[triListOffset + 1]];
        node->right->idxTriangle = m_triList[triListOffset + 1];
        return;
    }

    /* SAH splitting */

    enum class Axis { X, Y, Z };

    // Comparator for sorting triangles along an axis.
    auto triAxisGreater = [&](uint32_t t1, uint32_t t2, Axis axis) {
        if (axis == Axis::X && m_triAABBs[t1].min().x > m_triAABBs[t2].min().x)
            return true;
        else if (axis == Axis::Y && m_triAABBs[t1].min().y > m_triAABBs[t2].min().y)
            return true;
        else if (axis == Axis::Z && m_triAABBs[t1].min().z > m_triAABBs[t2].min().z)
            return true;
        return false;
        };

    float sahCost = std::numeric_limits<float>::max();
    Axis splitAxis = Axis::X;
    size_t splitPos = triListOffset + triCount / 2;

    // SAH: evaluate all split positions (between primitives) for each axis.
    // Cost: SA(L) * NL + SA(R) * NR
    for (const auto& axis : { Axis::X, Axis::Y, Axis::Z }) {
        std::sort(
            m_triList.begin() + triListOffset,
            m_triList.begin() + triListOffset + triCount,
            std::bind(triAxisGreater, std::placeholders::_1, std::placeholders::_2, axis)
        );

        // Prefix/suffix bounds to evaluate splits in O(n).
        std::vector<AABB> leftBounds(triCount);
        std::vector<AABB> rightBounds(triCount);

        leftBounds[0] = m_triAABBs[m_triList[triListOffset + 0]];
        for (size_t i = 1; i < triCount; i++) {
            leftBounds[i] = leftBounds[i - 1];
            leftBounds[i].merge(m_triAABBs[m_triList[triListOffset + i]]);
        }

        rightBounds[triCount - 1] = m_triAABBs[m_triList[triListOffset + triCount - 1]];
        for (size_t i = triCount - 1; i-- > 0;) {
            rightBounds[i] = rightBounds[i + 1];
            rightBounds[i].merge(m_triAABBs[m_triList[triListOffset + i]]);
        }

        // Split position is an index into m_triList.
        // Left: [offset, splitPos), Right: [splitPos, offset+count)
        for (size_t i = 1; i < triCount; i++) {
            float cost = leftBounds[i - 1].surfaceArea() * static_cast<float>(i);
            cost += rightBounds[i].surfaceArea() * static_cast<float>(triCount - i);
            if (cost < sahCost) {
                sahCost = cost;
                splitAxis = axis;
                splitPos = triListOffset + i;
            }
        }
    }

    if (splitPos <= triListOffset)
        splitPos = triListOffset + 1;
    else if (splitPos >= triListOffset + triCount)
        splitPos = triListOffset + triCount - 1;

    // Re-sort along the selected best axis so that splitPos corresponds to the final order.
    std::sort(
        m_triList.begin() + triListOffset,
        m_triList.begin() + triListOffset + triCount,
        std::bind(triAxisGreater, std::placeholders::_1, std::placeholders::_2, splitAxis)
    );

    /* Build children */
    node->left = std::make_unique<BvhNode>();
    buildRecursive(node->left.get(), triListOffset, splitPos - triListOffset);
    node->right = std::make_unique<BvhNode>();
    buildRecursive(node->right.get(), splitPos, triListOffset + triCount - splitPos);
}

std::vector<PathTracer::BufferBvhNode> PathTracer::BvhBufferizer::bufferize(BvhNode* root) {
    m_bufferData.clear();
    bufferizeRecursive(root);
    return m_bufferData;
}

void PathTracer::BvhBufferizer::bufferizeRecursive(BvhNode* node) {
    if (node == nullptr)
        return;
    BufferBvhNode bufferNode = {};
    bufferNode.idx = static_cast<uint32_t>(m_bufferData.size());
    bufferNode.aabbMin = Math::Vec4(node->aabb.min(), 0.0f);
    bufferNode.aabbMax = Math::Vec4(node->aabb.max(), 0.0f);
    if (node->left == nullptr && node->right == nullptr) {
        // Leaf node
        bufferNode.leafFlag = 1;
        bufferNode.idxTriangle = node->idxTriangle;
        m_bufferData.push_back(bufferNode);
    } else {
        // Internal node
        m_bufferData.push_back(bufferNode);
        bufferizeRecursive(node->left.get());
        m_bufferData[bufferNode.idx].rChildOffset =
            node->right != nullptr ?
            static_cast<uint32_t>(m_bufferData.size() - bufferNode.idx) :
            0; // 0 if no right child
        bufferizeRecursive(node->right.get());
    }
}
