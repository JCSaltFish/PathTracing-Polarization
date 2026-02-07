/**
 * @file PtScene.h
 * @brief Declaration of the PtScene class representing a scene containing multiple 3D models.
 */

#pragma once

#include "db/DbPub.h"
#include "utils/Math.h"

/**
 * @brief Represents a scene containing multiple 3D models.
 */
class PtScene {
    friend class DbTypeRegistry;
    /* OBJECT TYPE INFO */
private:
    static void serialize(DbSerializer& serializer, const PtScene& scene);
    static void deserialize(DbSerializer& serializer, PtScene& scene);
    static void migrate(int oldVersion, PtScene& scene);
public:
    static constexpr const char* TYPE_NAME = "PtScene";
    static constexpr int VERSION = 1;

    /* OBJECT FIELDS */
public:
    /**
     * @brief Camera settings for the scene.
    */
    struct Camera {
        Math::Vec3 position = Math::Vec3(0.0f, 0.0f, -10.0f); // Camera position in 3D space
        Math::Vec3 rotation = Math::Vec3(); // Camera rotation (Euler angles)
        float focusDist = 5.0f; // Focus distance
        float fStop = 32.0f; // F-stop value

        static constexpr float FOCAL = 0.05f; // Focal length (near clipping plane distance)
        static constexpr float FOV = 70.0f; // Field of view in degrees
    };
private:
    int m_traceDepth = 3; // Maximum recursion depth for path tracing
    int m_resX = 1024; // Image resolution X
    int m_resY = 768; // Image resolution Y

    Camera m_camera = {}; // Camera settings

    std::vector<DB::ID> m_models = {}; // List of model IDs in the scene

    /* OBJECT METHODS */
private:
    /**
     * @brief Get a const pointer to the PtScene object from a handle.
     * @param hScene Handle to the scene object.
     * @return Pointer to the PtScene object, or nullptr if invalid.
     */
    static const PtScene* view(const DbObjHandle& hScene);
public:
    /**
     * @brief Get the current trace depth.
     * @param hScene Handle to the scene object.
     * @return Current trace depth, or -1 if the handle is invalid.
     */
    static int getTraceDepth(const DbObjHandle& hScene);
    /**
     * @brief Set the trace depth.
     * @param hScene Handle to the scene object.
     * @param depth New trace depth value.
     * @return Result code indicating success or failure.
     */
    static DB::Result setTraceDepth(const DbObjHandle& hScene, int depth);
    /**
     * @brief Get the current resolution.
     * @param hScene Handle to the scene object.
     * @param[out] resX Reference to store the resolution X.
     * @param[out] resY Reference to store the resolution Y.
     */
    static void getResolution(const DbObjHandle& hScene, int& resX, int& resY);
    /**
     * @brief Set the resolution.
     * @param hScene Handle to the scene object.
     * @param resX New resolution X value.
     * @param resY New resolution Y value.
     * @return Result code indicating success or failure.
     */
    static DB::Result setResolution(const DbObjHandle& hScene, int resX, int resY);
    /**
     * @brief Get the camera settings.
     * @param hScene Handle to the scene object.
     * @return Camera settings, or default Camera if the handle is invalid.
     */
    static Camera getCamera(const DbObjHandle& hScene);
    /**
     * @brief Set the camera settings.
     * @param hScene Handle to the scene object.
     * @param camera New camera settings.
     * @return Result code indicating success or failure.
     */
    static DB::Result setCamera(const DbObjHandle& hScene, const Camera& camera);
    /**
     * @brief Get the list of models in the scene.
     * @param hScene Handle to the scene object.
     * @return Vector of handles to the models in the scene.
     */
    static std::vector<DbObjHandle> getModels(const DbObjHandle& hScene);
    /**
     * @brief Add a model to the scene.
     * @param hScene Handle to the scene object.
     * @param hModel Handle to the model object to add.
     * @return Result code indicating success or failure.
     */
    static DB::Result addModel(const DbObjHandle& hScene, const DbObjHandle& hModel);
    /**
     * @brief Remove a model from the scene.
     * @param hScene Handle to the scene object.
     * @param hModel Handle to the model object to remove.
     * @return Result code indicating success or failure.
     */
    static DB::Result delModel(const DbObjHandle& hScene, const DbObjHandle& hModel);
};
