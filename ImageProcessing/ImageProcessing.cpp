#define _CRT_SECURE_NO_WARNINGS

#include <iostream>
#include <string>
#include <vector>
#include <limits>
#include <thread>
#include <chrono>
#include <set>
#include "NaiveThreadPool.h"
#include "ImageData.h"
#include "atomic"
#include <sstream>
#include <map>
#include <algorithm>
#include <unordered_map>

#include "KdTree.h"

using TimerSignature = std::chrono::steady_clock::time_point;
using TimeStamp = std::chrono::duration<double>;

/** Atomic int for tracking processing progress */
static std::atomic<int> RowProgress;
static int AtomicSupport;

/** 
* Different strategies of calculating color similarity, each strategy can differ in speed drastically
* ColorPaletteComparison - calculates distance between colors - the lower the distance the more similar
* R_G_B_Sum - calculates sum of components (red, green, blue) - the closer the sums the more similar they are
*/
enum class EColorSimilarity : uint8_t
{
    Deprecated,
    R_G_B_Sum,
    Euclidean_Distance,
    Euclidean_Distance_KD_Tree
};

struct ColorRGB
{
    ColorRGB() : ColorRGB(0, 0, 0) {}

    ColorRGB(uint8_t r, uint8_t g, uint8_t b)
        : r(r), g(g), b(b)
    {

    }

    int GetComponentSum()
    {
        return r + g + b;
    }

    bool operator<(const ColorRGB& rhs) const
    {
        return r * 1000000 + g * 1000 + b < rhs.r * 1000000 + rhs.g * 1000 + rhs.b;
    }

    bool operator==(const ColorRGB& rhs)const
    {
        return r == rhs.r && g == rhs.g && b == rhs.b;
    }

    /** 
    * this hash function will cause sets and maps to grow very large, 
    * but at least there will be no collisions and its trivial to implement
    */
    struct HashFunction
    {
        size_t operator()(const ColorRGB& Color) const
        {
            return Color.r * 1000000 + Color.g * 1000 + Color.b;
        }
    };

    uint8_t r;
    uint8_t g;
    uint8_t b;
};

class SimpleScopedTimer 
{
public:
    SimpleScopedTimer(std::string Message)
        :Message(Message), StartTime(std::chrono::steady_clock::now())
    {
    }
    ~SimpleScopedTimer()
    {
        FinishTime = std::chrono::steady_clock::now();
        TimeStamp TotalTime = FinishTime - StartTime;
        std::cout << '\n' << Message << "Took: " << TotalTime.count() << " seconds to complete\n";
    }

    TimerSignature StartTime;
    TimerSignature FinishTime;

    std::string Message;
};

/** 
* Processes single row of pixels using R_G_B_Sum as strategy
* @see EColorSimilarity
*/
void taskProcessImage_RGB_Sum(ImageData& shapeImage, std::vector<ColorRGB>& ColorPalette, int row)
{
    int i = row;

    // pixels
    for (int j = 0; j < shapeImage.width; j++)
    {

        int pixelIndex = (i * shapeImage.width + j) * shapeImage.channel_num;
        constexpr std::size_t minSquaredDistance = std::numeric_limits<int>::max();

        // ptrs pointing to most similar color
        const unsigned char* minR = nullptr;
        const unsigned char* minG = nullptr;
        const unsigned char* minB = nullptr;
        //const unsigned char* minA = nullptr;

        // these are not necessary, but they make targetPixelColor assignment clearer
        unsigned char
            & r = shapeImage.data[pixelIndex],
            & g = shapeImage.data[pixelIndex + 1],
            & b = shapeImage.data[pixelIndex + 2];
            //& a = shapeImage.data[pixelIndex + 3];

        int targetPixelColor = r + g + b;

        // binary search through vector
        int left = 0;
        int right = ColorPalette.size() - 1;

        while (left <= right)
        {
            int middle = left + (right - left) / 2;
            if (ColorPalette[middle].GetComponentSum() == targetPixelColor)
            {
                // point to most similar color
                minR = &ColorPalette[middle].r;
                minG = &ColorPalette[middle].g;
                minB = &ColorPalette[middle].b;
                break;
            }
            else if (ColorPalette[middle].GetComponentSum() < targetPixelColor)
            {
                left = middle + 1;          
            }
            else 
            {
                right = middle - 1;
            }
        }

        // if the exact color was found in the palette array choose the most similar
        if (!minR)
        {
            int leftValue = std::abs(targetPixelColor - ColorPalette[left].GetComponentSum());
            int rightValue = std::abs(targetPixelColor - ColorPalette[left + 1].GetComponentSum());
            int mostSimilarColorIndex = leftValue > rightValue ? left + 1 : left;

            // point to most similar color
            minR = &ColorPalette[mostSimilarColorIndex].r;
            minG = &ColorPalette[mostSimilarColorIndex].g;
            minB = &ColorPalette[mostSimilarColorIndex].b;
        }

        // replace the pixel in-place
        shapeImage.data[pixelIndex] = *minR;
        shapeImage.data[pixelIndex + 1] = *minG;
        shapeImage.data[pixelIndex + 2] = *minB;
    }
    // update progress variable after processing single row
    ++RowProgress;
    AtomicSupport = RowProgress;
}

/** 
* replace every pixel in ShapeImage with appropriate/most similar color from ColorImage palette
* using map of pairs { ShapeImageColors - ColorImageColors }
*/
void ShapeImageToColorImagePalette(ImageData& shapeImage, 
    const std::unordered_map<ColorRGB, ColorRGB, ColorRGB::HashFunction>& ShapeToColorPaletteMap)
{
    // for every row
    for (int i = 0; i < shapeImage.height; i++)
    {
        // for every pixel
        for (int j = 0; j < shapeImage.width; j++)
        {
            // index of each pixel
            int pixelIndex = (i * shapeImage.width + j) * shapeImage.channel_num;

            int minSquaredDistance = std::numeric_limits<int>::max();

            // ptrs pointing to most similar color
            const unsigned char* minR = nullptr;
            const unsigned char* minG = nullptr;
            const unsigned char* minB = nullptr;

            // these are not necessary, but they make It "assignment" clearer
            unsigned char
                & r = shapeImage.data[pixelIndex],
                & g = shapeImage.data[pixelIndex + 1],
                & b = shapeImage.data[pixelIndex + 2];
                //& a = shapeImage.data[pixelIndex + 3];
            
            auto It = ShapeToColorPaletteMap.find(ColorRGB(r, g, b));

            // replace the pixel in-place
            shapeImage.data[pixelIndex] = It->second.r;
            shapeImage.data[pixelIndex + 1] = It->second.g;
            shapeImage.data[pixelIndex + 2] = It->second.b;
        }
        //// I don't think counter is neccesary here - the process seems to be quite fast compared to the rest of stuff
        //++RowProgress;
        //AtomicSupport = RowProgress;
    }
}


// calculate most similar color by comparing distance from each color from shape palette 
// to each color in color palette 
// complexity is O(n^2) - could be much better if I used octree or k-d tree
void taskCalculateColors_Euclidean_Distance(
    ImageData& ShapeImage,
    std::vector<ColorRGB>& ShapePalleteVec,
    std::set<ColorRGB>& ColorPalette,
    int i,
    std::unordered_map<ColorRGB, ColorRGB, ColorRGB::HashFunction>& Unordered_ShapeToColorPaletteMap,
    std::mutex& ColorPaletteMapMutex)
{
    // ptrs pointing to most similar color
    const unsigned char* minR = nullptr;
    const unsigned char* minG = nullptr;
    const unsigned char* minB = nullptr;
    /*const unsigned char* minA = nullptr;*/

    unsigned char
        & r = ShapePalleteVec[i].r,
        & g = ShapePalleteVec[i].g,
        & b = ShapePalleteVec[i].b;

    int minSquaredDistance = std::numeric_limits<int>::max();

    // calculate most similar color by comparing distance from each color from shape palette 
    // to each color in color palette O(n^2) - could be much better if I used octree or k-d tree
    for (const ColorRGB& colorPaletteColor : ColorPalette)
    {
        const unsigned char
            & r2 = colorPaletteColor.r,
            & g2 = colorPaletteColor.g,
            & b2 = colorPaletteColor.b;

        int squaredDistance = (r - r2) * (r - r2) + (g - g2) * (g - g2) + (b - b2) * (b - b2);

        if (squaredDistance < minSquaredDistance)
        {
            minSquaredDistance = squaredDistance;

            minR = &r2;
            minG = &g2;
            minB = &b2;
        }
    }
    ColorPaletteMapMutex.lock();
    Unordered_ShapeToColorPaletteMap.insert({ ShapePalleteVec[i], ColorRGB(*minR, *minG, *minB) });
    ColorPaletteMapMutex.unlock();
}

// stores unique color values of image in a set
std::set<ColorRGB> GetColorSetFromImage(ImageData image)
{
    std::set<ColorRGB> ColorSet;
    for (int i = 0; i < image.GetSizeInBytes(); i += image.channel_num)
    {
        ColorRGB color;
        color.r = image.data[i    ];
        color.g = image.data[i + 1];
        color.b = image.data[i + 2];
        ColorSet.insert(color);
    }
    return ColorSet;
}

void WriteImage(ImageData& ShapeImage, std::string ResultImageName)
{
    int stride = ShapeImage.channel_num * ShapeImage.width;
    //stride += (stride % 4) ? (4 - stride % 4) : 0;

    std::cout << "Loading Image...\n";

    int result = stbi_write_png(
        ResultImageName.c_str(),
        ShapeImage.width, ShapeImage.height,
        ShapeImage.channel_num, ShapeImage.data, stride);

    if (result)
    {
        std::cout << "Image writing Success!\n";
    }

    std::cout << "Job Finished!" << std::endl;
}

/**
* calculate sum of components for each color and pair it with the color itself:
* pair colors with their sums like this : { SumOfComponents - ColorRGB }
* @note that this will "cull" massive number of colors as this algorithm only allows maximum of 255 * 3 (765) colors
*
* then insert into vector for easier access
*/
void CalculateNewPaletteBasedOnComponentsSum(
    std::set<ColorRGB>& ColorPalette, 
    std::map<int, ColorRGB>& ColorPaletteSums, 
    std::vector<ColorRGB>& OutColorPaletteVector)
{
    for (const ColorRGB& color : ColorPalette)
    {
        ColorPaletteSums.insert({ color.r + color.g + color.b , color });
    }

    OutColorPaletteVector.reserve(ColorPaletteSums.size());
    for (const auto& pair : ColorPaletteSums)
    {
        OutColorPaletteVector.emplace_back(pair.second);
    }
}

void ConcurrentProcessImage_RGB_Sum(ImageData& ShapeImage, ImageData& ColorImage, 
    std::vector<ColorRGB>& ColorPaletteVector, int NumberOfWorkerThreads)
{
    ThreadPool Pool(NumberOfWorkerThreads);
    {
        for (int i = 0; i < ShapeImage.height; i++)
        {
            Pool.enqueue([&ShapeImage, &ColorImage, i, &ColorPaletteVector]()
                {
                    taskProcessImage_RGB_Sum(ShapeImage, ColorPaletteVector, i);
                });
        }

        while (Pool.busy())
        {
            std::stringstream msg;
            msg << "\rProcessed: " << AtomicSupport << " out of: " << ShapeImage.height << " rows";
            std::cout << msg.str();
        }
    }
}

void ConcurrentProcessImage_Euclidean_Distance(
    ImageData& ShapeImage, 
    std::vector<ColorRGB>& ShapePalleteVec,
    std::set<ColorRGB>& ColorPalette,
    std::unordered_map<ColorRGB, ColorRGB, ColorRGB::HashFunction>& Unordered_ShapeToColorPaletteMap,
    int NumberOfWorkerThreads)
{
    /**
    * mutex for locking acces to map
    * @note TODO I should come up with some way of accumulating results
    * and the adding them together instead of locking resource behind mutex
    */
    std::mutex ColorPaletteMapMutex;

    /**
    * TODO calculate for multiple colors per-task instead of calculating these for single color
    * this way there would be less queueing/locking/unlocking of thread pool most likely resulting in better performance
    */
    ThreadPool pool(NumberOfWorkerThreads);
    for (int i = 0; i < ShapePalleteVec.size(); i++)
    {
        pool.enqueue([&ShapeImage, &ShapePalleteVec, &ColorPalette, i, &Unordered_ShapeToColorPaletteMap, &ColorPaletteMapMutex]()
            {
                taskCalculateColors_Euclidean_Distance(
                    ShapeImage, ShapePalleteVec, ColorPalette,
                    i, Unordered_ShapeToColorPaletteMap, ColorPaletteMapMutex);
            });
    }
    while (pool.busy())
    {
        std::stringstream msg;
        msg << "\r" << Unordered_ShapeToColorPaletteMap.size();
        std::cout << msg.str();
    }
}

void PrepareKDTree(std::set<ColorRGB>& ColorPalette, KDTree<double>& tree) 
{
    std::vector<std::vector<double>> pixels;
    for (const ColorRGB& color : ColorPalette)
    {
        std::vector<double> planes{(double) color.r, (double) color.g, (double) color.b};
        pixels.push_back(planes);
    }
    tree.build(pixels);
}

void taskProcessImage_KD_Tree(ImageData& shapeImage, KDTree<double>& tree, int row)
{
    int i = row;

    // pixels
    for (int j = 0; j < shapeImage.width; j++)
    {
        int pixelIndex = (i * shapeImage.width + j) * shapeImage.channel_num;

        // these are not necessary, but they make targetPixelColor assignment clearer
        unsigned char
            & r = shapeImage.data[pixelIndex],
            & g = shapeImage.data[pixelIndex + 1],
            & b = shapeImage.data[pixelIndex + 2];
        //& a = shapeImage.data[pixelIndex + 3];

        std::vector<double> query{ (double)r, (double)g, (double)b };
        IndexedPoint<double> result = tree.nearestNeighbor(query);

        // replace the pixel in-place
        shapeImage.data[pixelIndex] = result.point[0];
        shapeImage.data[pixelIndex + 1] = result.point[1];
        shapeImage.data[pixelIndex + 2] = result.point[2];
    }
    // update progress variable after processing single row
    ++RowProgress;
    AtomicSupport = RowProgress;
}

void CalculateClosestColorNeighboursInKDTree(ImageData& ShapeImage, KDTree<double>& tree, int NumberOfWorkerThreads)
{
    ThreadPool Pool(NumberOfWorkerThreads);
    {
        for (int i = 0; i < ShapeImage.height; i++)
        {
            Pool.enqueue([&ShapeImage, i, &tree]()
                {
                    taskProcessImage_KD_Tree(ShapeImage, tree, i);
                });
        }

        while (Pool.busy())
        {
            std::stringstream msg;
            msg << "\rProcessed: " << AtomicSupport << " out of: " << ShapeImage.height << " rows";
            std::cout << msg.str();
        }
    }
}

void StartProcessingImage(ImageData& ShapeImage, ImageData& ColorImage, EColorSimilarity ImageProcessingStrategy, int NumberOfWorkerThreads)
{
    std::cout << "Hi, I process 2 images ShapeImage and ColorImage\n"
        << "As a result I redraw ShapeImage with color palette of ColorImage\n\n";

    /**
    * Load both images
    */
    std::cout << "Loading Images...\n";
    {
        SimpleScopedTimer bothImageTimer("Loading both images");

        {
            SimpleScopedTimer shapeImageTimer("Loading ShapeImage");
            ShapeImage.LoadImage();
        }
        {
            SimpleScopedTimer colorImageTimer("Loading ColorImage");
            ColorImage.LoadImage();
        }
    }

    std::cout << "Loading Images finished!\n\n";

    /** display pixel count for both images */
    std::cout << "Pixel count of shape image: " << ShapeImage.GetPixelCount() << " pixels\n";
    std::cout << "Pixel count of color image: " << ColorImage.GetPixelCount() << " pixels\n\n";

    /**
    * Calculate color palettes of both images
    */
    std::cout << "Loading color palettes...\n";
    std::set<ColorRGB> ShapePalette;
    std::set<ColorRGB> ColorPalette;
    {
        SimpleScopedTimer bothColorPalettesTimer("Loading both Palleters");
        {
            SimpleScopedTimer shapePaletteTimer("Loading ShapeImage Pallete");
            ShapePalette = GetColorSetFromImage(ShapeImage);
            std::cout << "ShapeImage has: " << ShapePalette.size() << " colors\n";
        }
        {
            SimpleScopedTimer colorPalleteTimer("Loading ColorImage Pallete");
            ColorPalette = GetColorSetFromImage(ColorImage);
            std::cout << "ColorImage has: " << ColorPalette.size() << " colors\n";
        }
    }

    std::cout << "Loading color palettes finished\n";
    std::string ResultImageName;
    {
        SimpleScopedTimer TotalProcessingTime("Image processing");
        switch (ImageProcessingStrategy)
        {
            case EColorSimilarity::R_G_B_Sum :
            {
                ResultImageName = "RGB_Sum_Result_Stbi.png";
                std::map<int, ColorRGB> ColorPaletteSums;
                std::vector<ColorRGB> ColorPaletteVector;
                CalculateNewPaletteBasedOnComponentsSum(ColorPalette, ColorPaletteSums, ColorPaletteVector);
                ConcurrentProcessImage_RGB_Sum(ShapeImage, ColorImage, ColorPaletteVector, NumberOfWorkerThreads);
                break;
            }

            case EColorSimilarity::Euclidean_Distance :
            {
                ResultImageName = "Euclidean_Distance_Result_Stbi.png";
                std::unordered_map<ColorRGB, ColorRGB, ColorRGB::HashFunction> Unordered_ShapeToColorPaletteMap;
                std::vector<ColorRGB> ShapePalleteVec(ShapePalette.begin(), ShapePalette.end());

                // create a map of pairs { colorRGB(ShapeImage) - Most Similar colorRGB(ColorImage) }
                ConcurrentProcessImage_Euclidean_Distance(ShapeImage, ShapePalleteVec, 
                    ColorPalette, Unordered_ShapeToColorPaletteMap, NumberOfWorkerThreads);

                // apply map to ShapeImage
                ShapeImageToColorImagePalette(ShapeImage, Unordered_ShapeToColorPaletteMap);
                break;
            }

            case EColorSimilarity::Euclidean_Distance_KD_Tree : 
            {
                ResultImageName = "Euclidean_Distance_KD_Tree_Result_Stbi.png";
                KDTree<double> tree;

                PrepareKDTree(ColorPalette, tree);
                CalculateClosestColorNeighboursInKDTree(ShapeImage, tree, NumberOfWorkerThreads);
                break;
            }
        }
    }

    {
        SimpleScopedTimer TotalProcessingTime("total writing result image");
        WriteImage(ShapeImage, ResultImageName);
    }
}

int main()
{
    ImageData ShapeImage("ShapeImage.png");
    ImageData ColorImage("ColorImage.png");
    EColorSimilarity ImageProcessingStrategy = EColorSimilarity::Euclidean_Distance_KD_Tree;
    int NumberOfWorkerThreads = 31;

    {
        SimpleScopedTimer TotalProcessingTime("total time of program running");
        StartProcessingImage(ShapeImage, ColorImage, ImageProcessingStrategy, NumberOfWorkerThreads);
    }
}