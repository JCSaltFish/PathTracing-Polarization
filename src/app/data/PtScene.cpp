/**
 * @file PtScene.cpp
 * @brief Implementation of the PtScene class representing a scene containing multiple 3D models.
 */

#include "app/AppDataManager.h"

void PtScene::serialize(DbSerializer& serializer, const PtScene& scene) {
    serializer.serialize(scene.m_traceDepth);
    serializer.serialize(scene.m_resX);
    serializer.serialize(scene.m_resY);

    serializer.serialize(scene.m_camera.position.x);
    serializer.serialize(scene.m_camera.position.y);
    serializer.serialize(scene.m_camera.position.z);

    serializer.serialize(scene.m_camera.rotation.x);
    serializer.serialize(scene.m_camera.rotation.y);
    serializer.serialize(scene.m_camera.rotation.z);

    serializer.serialize(scene.m_camera.focusDist);
    serializer.serialize(scene.m_camera.fStop);

    serializer.serialize(scene.m_models);
}

void PtScene::deserialize(DbSerializer& serializer, PtScene& scene) {
    serializer.deserialize(scene.m_traceDepth);
    serializer.deserialize(scene.m_resX);
    serializer.deserialize(scene.m_resY);

    serializer.deserialize(scene.m_camera.position.x);
    serializer.deserialize(scene.m_camera.position.y);
    serializer.deserialize(scene.m_camera.position.z);

    serializer.deserialize(scene.m_camera.rotation.x);
    serializer.deserialize(scene.m_camera.rotation.y);
    serializer.deserialize(scene.m_camera.rotation.z);

    serializer.deserialize(scene.m_camera.focusDist);
    serializer.deserialize(scene.m_camera.fStop);

    serializer.deserialize(scene.m_models);
}

void PtScene::migrate(int oldVersion, PtScene& scene) {}

const PtScene* PtScene::view(const DbObjHandle& hScene) {
    if (!hScene.isValid() || hScene.getType() != PtScene::TYPE_NAME)
        return nullptr;
    return hScene.getDB()->objGet<PtScene>(hScene);
}

int PtScene::getTraceDepth(const DbObjHandle& hScene) {
    const PtScene* scene = view(hScene);
    if (!scene)
        return 0;
    return scene->m_traceDepth;
}

DB::Result PtScene::setTraceDepth(const DbObjHandle& hScene, int depth) {
    const PtScene* scene = view(hScene);
    if (!scene)
        return DB::Result::INVALID_HANDLE;
    if (depth < 0)
        return DB::Result::FAILURE;
    PtScene newScene = *scene;
    newScene.m_traceDepth = depth;
    return hScene.getDB()->objModify(hScene, newScene);
}

void PtScene::getResolution(const DbObjHandle& hScene, int& resX, int& resY) {
    const PtScene* scene = view(hScene);
    if (!scene) {
        resX = 0;
        resY = 0;
        return;
    }
    resX = scene->m_resX;
    resY = scene->m_resY;
}

DB::Result PtScene::setResolution(const DbObjHandle& hScene, int resX, int resY) {
    const PtScene* scene = view(hScene);
    if (!scene)
        return DB::Result::INVALID_HANDLE;
    if (resX <= 0 || resY <= 0)
        return DB::Result::FAILURE;
    PtScene newScene = *scene;
    newScene.m_resX = resX;
    newScene.m_resY = resY;
    return hScene.getDB()->objModify(hScene, newScene);
}

PtScene::Camera PtScene::getCamera(const DbObjHandle& hScene) {
    const PtScene* scene = view(hScene);
    if (!scene)
        return {};
    return scene->m_camera;
}

DB::Result PtScene::setCamera(const DbObjHandle& hScene, const Camera& camera) {
    const PtScene* scene = view(hScene);
    if (!scene)
        return DB::Result::INVALID_HANDLE;
    PtScene newScene = *scene;
    newScene.m_camera = camera;
    return hScene.getDB()->objModify(hScene, newScene);
}

std::vector<DbObjHandle> PtScene::getModels(const DbObjHandle& hScene) {
    std::vector<DbObjHandle> models;
    const PtScene* scene = view(hScene);
    if (!scene)
        return models;
    DB* db = hScene.getDB();
    for (DB::ID modelId : scene->m_models) {
        DbObjHandle hModel(db, modelId);
        if (hModel.isValid() && hModel.getType() == PtModel::TYPE_NAME)
            models.push_back(hModel);
    }
    return models;
}

DB::Result PtScene::addModel(const DbObjHandle& hScene, const DbObjHandle& hModel) {
    const PtScene* scene = view(hScene);
    if (!scene)
        return DB::Result::INVALID_HANDLE;
    if (!hModel.isValid() || hModel.getType() != PtModel::TYPE_NAME)
        return DB::Result::INVALID_HANDLE;

    if (std::binary_search(scene->m_models.begin(), scene->m_models.end(), hModel.getID()))
        return DB::Result::SUCCESS; // Already in the scene

    PtScene newScene = *scene;
    newScene.m_models.push_back(hModel.getID());
    return hScene.getDB()->objModify(hScene, newScene);
}

DB::Result PtScene::delModel(const DbObjHandle& hScene, const DbObjHandle& hModel) {
    const PtScene* scene = view(hScene);
    if (!scene)
        return DB::Result::INVALID_HANDLE;
    if (!hModel.isValid() || hModel.getType() != PtModel::TYPE_NAME)
        return DB::Result::INVALID_HANDLE;

    if (!std::binary_search(scene->m_models.begin(), scene->m_models.end(), hModel.getID()))
        return DB::Result::SUCCESS; // Not in the scene

    PtScene newScene = *scene;
    newScene.m_models.erase
    (
        std::remove(newScene.m_models.begin(), newScene.m_models.end(), hModel.getID()),
        newScene.m_models.end()
    );
    DB::Result result = hScene.getDB()->objModify(hScene, newScene);
    if (result != DB::Result::SUCCESS)
        return result;

    DB* db = hScene.getDB();
    // Also delete the model and its meshes and materials
    for (auto& hMesh : PtModel::getMeshes(hModel)) {
        DbObjHandle hMaterial = PtMesh::getMaterial(hMesh);
        if (hMaterial.isValid()) {
            result = db->objDelete<PtMaterial>(hMaterial);
            if (result != DB::Result::SUCCESS)
                return result;
        }
        result = db->objDelete<PtMesh>(hMesh);
        if (result != DB::Result::SUCCESS)
            return result;
    }
    return db->objDelete<PtModel>(hModel);
}
