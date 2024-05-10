#include <iostream>
#include <fstream>
#include <unistd.h>
#include <math.h>
#include <vector>
#include "Obj.h"

#define X 0
#define Y 1
#define Z 2
#define MUR

#define VECTOR std::vector<float>
#define VECTOR_2D std::vector<std::vector<float>>
#define VECTOR_3D std::vector<std::vector<std::vector<float>>>
#define VECTOR_4D std::vector<std::vector<std::vector<std::vector<float>>>>

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

// Для граничных условий мура
const float Mur_Tx = c * dt / dx;
const float Mur_Ty = c * dt / dy;
const float Mur_Tz = c * dt / dz;

struct Mesh3D
{
    int len_x, len_y, len_z;
};

class Field
{
public:
    Field() {}
    Field(Mesh3D &mesh, Obj &obj, float value = 0)
    {
        len_x = mesh.len_x;
        len_y = mesh.len_y;
        len_z = mesh.len_z;
        obj_x = obj.get_obj_x();
        obj_y = obj.get_obj_y();
        obj_z = obj.get_obj_z();
        _obj_len_x = obj.get__obj_len_x();
        _obj_len_y = obj.get__obj_len_y();
        _obj_len_z = obj.get__obj_len_z();
        this->_value = value;
        nodes = VECTOR_4D(len_x,
                          VECTOR_3D(len_y,
                                    VECTOR_2D(len_z, VECTOR(3, _value))));
        nodes_old = VECTOR_4D(len_x,
                              VECTOR_3D(len_y,
                                        VECTOR_2D(len_z, VECTOR(3, _value))));
        nodesE = VECTOR_3D(len_x, VECTOR_2D(len_y, VECTOR(3, _value)));
    }

    float gt(float x, float y, float z, int comp);
    void setNode(float x, float y, float z, int comp, float value);
    float gtO(float x, float y, float z, int comp);

    void updateE(Field &H, Field &J, Field &c1, Field &c2);
    void updateH(Field &E, Field &c);

    void saveExToFileZ(const std::string &filename, int comp);
    void saveExToFileX(const std::string &filename, int comp);
    void saveExToFileY(const std::string &filename, int comp);
    void saveToFile(const std::string &filename);

    void saveGraphZ(const std::string &filename, int comp);
    void saveRangeY();

    void copyNodes();

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
    int obj_x = 0;
    int obj_y = 0;
    int obj_z = 0;
    int _obj_len_x = 0;
    int _obj_len_y = 0;
    int _obj_len_z = 0;
    float _value = 0;
    float min_x = 0;
    float min_y = 0;
    float min_z = 0;
    float max_x = 0;
    float max_y = 0;
    float max_z = 0;
    VECTOR_4D nodes;
    VECTOR_4D nodes_old;
    VECTOR_3D nodesE;
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
    int obj_x = 75;
    int obj_y = 75;
    int _obj_len_x = 50;
    int _obj_len_y = 50;
    Field E;
    Field H;
    Field Sigm; // 2 -14
    Field J;
    Field Eps; // 1.00057
    Field Mu;  // 0.999991

    FDTD(Mesh3D &mesh, Obj &obj)
    {
        E = Field(mesh, obj);
        H = Field(mesh, obj);
        J = Field(mesh, obj);
        Sigm = Field(mesh, obj, 2.0e-5); // 2 -14
        Eps = Field(mesh, obj, 1.00057);  // 1.00057
        Mu = Field(mesh, obj, 0.999991);  // 0.999991
        c = Field(mesh, obj);
        c1 = Field(mesh, obj);
        c2 = Field(mesh, obj);
        len_x = mesh.len_x;
        len_y = mesh.len_y;
        len_z = mesh.len_z;
        initCoeffi();
    }
    void update(float time);
    void initCoeffi();
    void Mur();
    void GetBorderValues();
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
//                 float amplitude = (amplitudes.gt(x, y, z, X) + amplitudes.gt(x, y, z, Y) + amplitudes.gt(x, y, z, Z)) / 3;

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