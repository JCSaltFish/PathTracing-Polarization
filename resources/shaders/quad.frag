/**
 * @file quad.frag
 * @brief Fragment shader for rendering a textured quad.
 */

#version 450

layout(location = 0) in vec2 v_texCoord; // Texture coordinate

layout(binding = 0) uniform sampler2D u_texture; // Texture sampler

/**
 * @brief Uniform buffer for shader parameters.
 */
layout(binding = 1) uniform Params {
    int channel; // Channel selection (0 for red, 1 for green, 2 for blue)
} u_params; // Shader parameters

layout(location = 0) out vec4 o_fragColor; // Final fragment color

void main() {
    vec3 color = texture(u_texture, v_texCoord).rgb;
    o_fragColor = vec4(0.0, 0.0, 0.0, 1.0);
    if (u_params.channel == 0) {
        color *= 5.0;
        o_fragColor = vec4(color.rrr, 1.0);
    } else {
        color = 1.0 - (color + 0.025) * 20.0;
        if (u_params.channel == 1)
            o_fragColor = vec4(color.ggg, 1.0);
        else if (u_params.channel == 2)
            o_fragColor = vec4(color.bbb, 1.0);
    }
}
