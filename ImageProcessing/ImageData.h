#pragma once
#include <string>
// stb includes
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
// end stb includes
struct ImageData
{
    ImageData(std::string newImagePath)
        : imagePath(newImagePath)
    {

    }

    bool LoadImage()
    {
        data = stbi_load(imagePath.c_str(), &width, &height, &channel_num, 0);
        return true;
    }

    int GetPixelCount()
    {
        return width * height;
    }

    int GetSizeInBytes()
    {
        return GetPixelCount() * channel_num;
    }

    std::string imagePath;
    int width;
    int height;
    int channel_num;
    unsigned char* data;
};