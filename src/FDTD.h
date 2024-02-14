#include <iostream>
#include <fstream>
#include <unistd.h>
#include <math.h>
#include <vector>

#define X 0
#define Y 1
#define Z 2

// Константы и параметры для моделирования
const float c = 2.997925010e8;       // Скорость света (м/с)
const float epsilon0 = 8.854e-12;    // Вакуумная диэлектрическая проницаемость (Ф/м)
const float mu0 = 4 * M_PI * 1.0e-7; // Вакуумная магнитная проницаемость (Гн/м)
const float dt = 1.0e-12;            // Временной шаг (с)
const float dstep = 0.3e-2;
const float dx = dstep; // Пространственный шаг по x (м)
const float dy = dstep; // Пространственный шаг по y (м)
const float dz = dstep; // Пространственный шаг по z (м)
const float sigma = 0;  // Коэффициент проводимости для потерь

struct Mesh3D
{
    int len_x, len_y, len_z;
};

class Field
{
public:
    Field(){}
    Field(Mesh3D &mesh, float value = 0)
    {
        len_x = mesh.len_x;
        len_y = mesh.len_y;
        len_z = mesh.len_z;
        this->_value = value;
        nodes = std::vector<std::vector<std::vector<std::vector<float>>>>(len_x,
                                                                          std::vector<std::vector<std::vector<float>>>(len_y,
                                                                                                                       std::vector<std::vector<float>>(len_z, std::vector<float>(3, _value))));
        nodesE = std::vector<std::vector<std::vector<float>>>(len_x, std::vector<std::vector<float>>(len_y, std::vector<float>(3, _value)));
    }

    float getNode(float x, float y, float z, int comp);

    float getNodeE(float x, float y, int comp);

    void setNode(float x, float y, float z, int comp, float value);

    void updateE(Field &H, Field &J, Field &c1 , Field &c2);
    void updateE();

    void updateH(Field &E, Field &c);

    void saveExToFileZ(const std::string &filename, int comp);
    void saveExToFileX(const std::string &filename, int comp);
    void saveExToFileY(const std::string &filename, int comp);

    void saveGraphZ(const std::string &filename, int comp);
    void saveRangeY();

    size_t getSizeX()
    {
        return len_x;
    }
    size_t getSizeY()
    {
        return len_y;
    }
    size_t getSizeZ()
    {
        return len_z;
    }
    float getMinX() const
    {
        return min_x;
    }

    float getMinY() const
    {
        return min_y;
    }

    float getMinZ() const
    {
        return min_z;
    }

    float getMaxX() const
    {
        return max_x;
    }

    float getMaxY() const
    {
        return max_y;
    }

    float getMaxZ() const
    {
        return max_z;
    }

private:
    int len_x = 0;
    int len_y = 0;
    int len_z = 0;
    float _value = 0;
    float min_x = 0;
    float min_y = 0;
    float min_z = 0;
    float max_x = 0;
    float max_y = 0;
    float max_z = 0;
    std::vector<std::vector<std::vector<std::vector<float>>>> nodes;
    std::vector<std::vector<std::vector<float>>> nodesE;
};

class FDTD
{
private:
    Field c;
    Field c1;
    Field c2;
    float _time = 0;
    int len_x = 0;
    int len_y = 0;
    int len_z = 0;

public:
    Field E;
    Field H;
    Field Sigm; // 2 -14
    Field J;
    Field Eps; // 1.00057
    Field Mu;  // 0.999991

    FDTD(Mesh3D &mesh)
    {
        E = Field(mesh);
        H = Field(mesh);
        J = Field(mesh);
        Sigm = Field(mesh, 2.0e-14); // 2 -14
        Eps = Field(mesh, 1.00057);  // 1.00057
        Mu = Field(mesh, 0.999991);  // 0.999991
        c = Field(mesh);
        c1 = Field(mesh);
        c2 = Field(mesh);
        len_x = mesh.len_x;
        len_y = mesh.len_y;
        len_z = mesh.len_z;
        initCoeffi();
    }
    void update(float time);
    void initCoeffi();
};





// void drawWaves(sf::RenderWindow &window, Field &amplitudes)
// {
//     // Размер окна
//     const int width = window.getSize().x;
//     const int height = window.getSize().y;
//     size_t x_size = amplitudes.getSizeX();
//     size_t y_size = amplitudes.getSizeY();
//     size_t z_size = amplitudes.getSizeZ();

//     // Создаем изображение
//     sf::Image image;
//     image.create(width, height, sf::Color::White);

//     // Проходим по каждой точке изображения и устанавливаем цвет в зависимости от амплитуды
//     for (int x = 0; x < x_size; ++x)
//     {
//         for (int y = 0; y < y_size; ++y)
//         {
//             for (int z = 0; z < z_size; ++z)
//             {
//                 float amplitude = (amplitudes.getNode(x, y, z, X) + amplitudes.getNode(x, y, z, Y) + amplitudes.getNode(x, y, z, Z)) / 3;

//                 sf::Color color(static_cast<sf::Uint8>(amplitude * 255), 0, 0);
//                 image.setPixel(x, y, color);
//             }
//         }
//     }

//     // Создаем текстуру изображения
//     sf::Texture texture;
//     texture.loadFromImage(image);

//     // Создаем спрайт с текстурой
//     sf::Sprite sprite(texture);

//     // Отрисовываем спрайт
//     window.draw(sprite);
// }