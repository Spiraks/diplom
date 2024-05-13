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
#define PML

#define VECTOR std::vector<float>
#define VECTOR_2D std::vector<VECTOR>
#define VECTOR_3D std::vector<VECTOR_2D>
#define VECTOR_4D std::vector<VECTOR_3D>

#define V_PML std::vector<PMLData>
#define V_PML_2D std::vector<V_PML>
#define V_PML_3D std::vector<V_PML_2D>

// Константы и параметры для моделирования
const float c = 2.997925010e8;   //Скорость света (м/с)
const float eps0 = 8.854e-12;    //Вакуумная диэлектрическая проницаемость (Ф/м)
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

struct PMLData {
    float Hyx, Hyz, Hzx, Hzy, Hxy, Hxz;
    float Exy, Exz, Eyx, Eyz, Ezx, Ezy;
};

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
        nodes = VECTOR_4D(len_x, VECTOR_3D(len_y, VECTOR_2D(len_z, VECTOR(3, _value))));
        nodes_1 = VECTOR_4D(len_x, VECTOR_3D(len_y, VECTOR_2D(len_z, VECTOR(3, _value))));
        nodesE = VECTOR_3D(len_x, VECTOR_2D(len_y, VECTOR(3, _value)));
    }

    float gt(size_t x, size_t y, size_t z, size_t comp);
    void setNode(size_t x, size_t y, size_t z, size_t comp, float value);
    float gtO(size_t x, size_t y, size_t z, size_t comp);

    void updateE(Field &H, Field &J, Field &c1, Field &c2);
    void updateH(Field &E, Field &c0);

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
    //for mur
    VECTOR_4D nodes_1;
    VECTOR_3D nodesE;
};

class FDTD
{
private:
    Field c0;
    Field c1;
    Field c2;
    float _time = 0;
    int len_x = 0;
    int len_y = 0;
    int len_z = 0;

    //for pml

        // по граням
    V_PML_3D pmlXN;
    V_PML_3D pmlX ;
    V_PML_3D pmlYN;
    V_PML_3D pmlY ;
    V_PML_3D pmlZN;
    V_PML_3D pmlZ ;
    // по ребрам
    // по углам

public:
    int obj_x = 75;
    int obj_y = 75;
    int _obj_len_x = 50;
    int _obj_len_y = 50;

    size_t s = 1;
    size_t s1 = s + 1;
    size_t end = 2;
    size_t x = len_x - end;
    size_t y = len_y - end;
    size_t z = len_z - end;
    size_t x1 = len_x - end - 1;
    size_t y1 = len_y - end - 1;
    size_t z1 = len_z - end - 1;
    //for pml
    float pmlBHX = 1;
    float pmlBHZ = 1;
    float pmlBHY = 1;
    int plmLayerNumber = 10;
    size_t plmLayerNumberMinus1 = plmLayerNumber - 1;
    size_t plmLayerNumberMinus2 = plmLayerNumber - 2;
    std::vector<float> pmlSigmaStarH;
    std::vector<float> pmlSigmaStarE;
    std::vector<float> pmlExpSigmaStarX; 

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
        c0 = Field(mesh, obj);
        c1 = Field(mesh, obj);
        c2 = Field(mesh, obj);
        len_x = mesh.len_x;
        len_y = mesh.len_y;
        len_z = mesh.len_z;
        pmlSigmaStarH = VECTOR(plmLayerNumber);
        pmlSigmaStarE = VECTOR(plmLayerNumber);
        pmlExpSigmaStarX = VECTOR(plmLayerNumber,1);



        pmlXN = V_PML_3D(plmLayerNumber,V_PML_2D(len_y, V_PML(len_z)));
        pmlX  = V_PML_3D(plmLayerNumber,V_PML_2D(len_y, V_PML(len_z)));
        pmlYN = V_PML_3D(len_x,V_PML_2D(plmLayerNumber, V_PML(len_z)));
        pmlY  = V_PML_3D(len_x,V_PML_2D(plmLayerNumber, V_PML(len_z)));
        pmlZN = V_PML_3D(len_x,V_PML_2D(len_y, V_PML(plmLayerNumber)));
        pmlZ  = V_PML_3D(len_x,V_PML_2D(len_y, V_PML(plmLayerNumber)));
        
        initCoeffi();
    }
    void update(float time);
    void initCoeffi();

    void Mur();
    void GetBorderValuesMur();

    void Pml();
    void PML_E();
    void PML_H();
    void InitPML();
    float GetSigma(float D, float dStep, float pmlG, float pmlR,int N);

    //подсчет граней
    void PML_UpdateEXN();
    void PML_UpdateEX();
    void PML_UpdateEYN();
    void PML_UpdateEY();
    void PML_UpdateEZN();
    void PML_UpdateEZ();
    //подсчет ребер
    // void PML_UpdateEXNYN();
    // void PML_UpdateEXYN();
    // void PML_UpdateEXNY();
    // void PML_UpdateEXY();

    // void PML_UpdateEZNYN();
    // void PML_UpdateEZYN();
    // void PML_UpdateEZNY();
    // void PML_UpdateEZY();

    // void PML_UpdateEZNXN();
    // void PML_UpdateEZXN();
    // void PML_UpdateEZNX();
    // void PML_UpdateEZX();

    // // подсчет углов
    // void PML_UpdateEXYZ();
    // void PML_UpdateEXNYZ();
    // void PML_UpdateEXYZN();
    // void PML_UpdateEXYNZ();

    // void PML_UpdateEXNYZN();
    // void PML_UpdateEXYNZN();
    // void PML_UpdateEXNYNZ();
    // void PML_UpdateEXNYNZN();

    //подсчет граней
    void PML_UpdateHXN();
    void PML_UpdateHX();
    void PML_UpdateHYN();
    void PML_UpdateHY();
    void PML_UpdateHZN();
    void PML_UpdateHZ();
    //подсчет ребер
    // void PML_UpdateHXNYN();
    // void PML_UpdateHXYN();
    // void PML_UpdateHXNY();
    // void PML_UpdateHXY();

    // void PML_UpdateHZNYN();
    // void PML_UpdateHZYN();
    // void PML_UpdateHZNY();
    // void PML_UpdateHZY();

    // void PML_UpdateHZNXN();
    // void PML_UpdateHZXN();
    // void PML_UpdateHZNX();
    // void PML_UpdateHZX();

    // // подсчет углов
    // void PML_UpdateHXYZ();
    // void PML_UpdateHXNYZ();
    // void PML_UpdateHXYZN();
    // void PML_UpdateHXYNZ();

    // void PML_UpdateHXNYZN();
    // void PML_UpdateHXYNZN();
    // void PML_UpdateHXNYNZ();
    // void PML_UpdateHXNYNZN();
};