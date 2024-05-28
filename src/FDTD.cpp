
#include "FDTD.h"
#include <cmath>


void Field::copyNodes()
{
    _X_1 = _X;
    _Y_1 = _Y;
    _Z_1 = _Z;
};
void Field::updateE(Field &H, Field &J, VECTOR_3D &AE, VECTOR_3D &BE)
{

    for (size_t z = 1; z < len_z - 1; z++)
    {
        for (size_t y = 1; y < len_y - 1; y++)
        {
            for (size_t x = 1; x < len_x - 1; x++)
            {
                _X[x][y][z]=_X[x][y][z] * AE[x][y][z]+((H._Z[x][y+1][z]-H._Z[x][y-1][z])/dy-(H._Y[x][y][z+1]-H._Y[x][y][z-1])/dz)*BE[x][y][z] - J._X[x][y][z];

                _Y[x][y][z]=_Y[x][y][z] * AE[x][y][z]+((H._X[x][y][z+1]-H._X[x][y][z-1])/dz-(H._Z[x+1][y][z]-H._Z[x-1][y][z])/dx)*BE[x][y][z]  - J._Y[x][y][z];

                _Z[x][y][z]=_Z[x][y][z] * AE[x][y][z]+((H._Y[x+1][y][z]-H._Y[x-1][y][z])/dx-(H._X[x][y+1][z]-H._X[x][y-1][z])/dy)*BE[x][y][z] - J._Z[x][y][z];
            }
        }
    }
}
int t = 0;
void Field::updateH(Field &E, VECTOR_3D &AH)
{
    for (size_t z = 1; z < len_z - 1; z++)
    {
        for (size_t y = 1; y < len_y - 1; y++)
        {
            for (size_t x = 1; x < len_x - 1; x++)
            {
                _X[x][y][z]=_X[x][y][z]-((E._Z[x][y+1][z]-E._Z[x][y-1][z])/dy-(E._Y[x][y][z+1]-E._Y[x][y][z-1])/dz)*AH[x][y][z];

                _Y[x][y][z] =  _Y[x][y][z]-((E._X[x][y][z+1]-E._X[x][y][z-1])/dz-(E._Z[x+1][y][z]-E._Z[x-1][y][z])/dx)*AH[x][y][z];

                _Z[x][y][z]=_Z[x][y][z]-((E._Y[x+1][y][z]-E._Y[x-1][y][z])/dx-(E._X[x][y+1][z]-E._X[x][y-1][z])/dy)*AH[x][y][z];
            }
        }
    }
}
static const int slice = 2;
#ifdef FLEX
void Field::saveExToFileZX(const std::string &filename)
{
    std::ofstream file(filename);

    if (!file.is_open())
    {
        std::cerr << "Ошибка открытия файла для записи: " << filename << std::endl;
        return;
    }
    int z = len_z / slice;
    for (size_t y = 1; y < len_y-1; y++)
    {
        for (size_t x = 1; x < len_x-1; x++)

        {
            float val = _X[x][y][z];
            
            if (val == 0)
            {
                val = _X[x-1]   [y]     [z] +
                      _X[x]     [y-1]   [z] +
                      _X[x+1]   [y]     [z] +
                      _X[x]     [y+1]   [z] +
                      _X[x+1]   [y+1]   [z] +
                      _X[x+1]   [y-1]   [z] +
                      _X[x-1]   [y+1]   [z] +
                      _X[x-1]   [y-1]   [z];
                val = val / 8;
            }
            
            file << val << " ";
        }
        file << std::endl;
    }
    file.close();
    std::cout << "Данные EX сохранены в файле: " << filename << std::endl;
}

void Field::saveExToFileZY(const std::string &filename)
{
    std::ofstream file(filename);

    if (!file.is_open())
    {
        std::cerr << "Ошибка открытия файла для записи: " << filename << std::endl;
        return;
    }
    int z = len_z / slice;
    for (size_t y = 1; y < len_y-1; y++)
    {
        for (size_t x = 1; x < len_x-1; x++)

        {
            float val = _Y[x][y][z];
            
            if (val == 0)
            {
                val = _Y[x - 1] [y]     [z] +
                      _Y[x]     [y - 1] [z] +
                      _Y[x + 1] [y]     [z] +
                      _Y[x]     [y + 1] [z] +
                      _Y[x + 1] [y + 1] [z] +
                      _Y[x + 1] [y - 1] [z] +
                      _Y[x - 1] [y + 1] [z] +
                      _Y[x - 1] [y - 1] [z];
                val = val / 8;
            }
            
            file << val << " ";
        }
        file << std::endl;
    }

    file.close();
    std::cout << "Данные EY сохранены в файле: " << filename << std::endl;
}

void Field::saveExToFileYX(const std::string &filename)
{
    std::ofstream file(filename);

    if (!file.is_open())
    {
        std::cerr << "Ошибка открытия файла для записи: " << filename << std::endl;
        return;
    }
    int y = len_y / slice;
    for (size_t x = 1; x < len_x-1; x++)
    {
    for (size_t z = 1; z < len_z-1; z++)
        {
            float val = _X[x][y][z];
            
            if (val == 0)
            {
                val = _X[x-1]   [y] [z]     +
                      _X[x]     [y] [z-1]   +
                      _X[x+1]   [y] [z]     +
                      _X[x]     [y] [z+1]   +
                      _X[x+1]   [y] [z+1]   +
                      _X[x+1]   [y] [z-1]   +
                      _X[x-1]   [y] [z+1]   +
                      _X[x-1]   [y] [z-1]   ;
                val = val / 8;
            }
            
            file << val << " ";
        }
        file << std::endl;
    }

    file.close();
    std::cout << "Данные EX сохранены в файле: " << filename << std::endl;
}

void Field::saveExToFileYZ(const std::string &filename)
{
    std::ofstream file(filename);

    if (!file.is_open())
    {
        std::cerr << "Ошибка открытия файла для записи: " << filename << std::endl;
        return;
    }
    int y = len_y / slice;
    for (size_t x = 1; x < len_x-1; x++)

    {
    for (size_t z = 1; z < len_z-1; z++)
        {
            float val = _Z[x][y][z];
            
            if (val == 0)
            {
                val =_Z[x-1]   [y] [z]     +
                    _Z[x]     [y] [z-1]   +
                    _Z[x+1]   [y] [z]     +
                    _Z[x]     [y] [z+1]   +
                    _Z[x+1]   [y] [z+1]   +
                    _Z[x+1]   [y] [z-1]   +
                    _Z[x-1]   [y] [z+1]   +
                    _Z[x-1]   [y] [z-1]   ;
                val = val / 8;
            }
            
            file << val << " ";
        }
        file << std::endl;
    }

    file.close();
    std::cout << "Данные EZ сохранены в файле: " << filename << std::endl;
}

void Field::saveExToFileXZ(const std::string &filename)
{
    std::ofstream file(filename);

    if (!file.is_open())
    {
        std::cerr << "Ошибка открытия файла для записи: " << filename << std::endl;
        return;
    }
    int x = len_x / slice;
    for (size_t y = 1; y < len_y-1; y++)
    {
    for (size_t z = 1; z < len_z-1; z++)
        {
            float val = _Z[x][y][z];
            
            if (val == 0)
            {
                val = _Z[x][y][z-1] +
                      _Z[x][y-1][ z] +
                      _Z[x][ y][ z + 1] +
                      _Z[x][ y + 1][ z] +
                      _Z[x][ y + 1][ z + 1] +
                      _Z[x][ y - 1][ z + 1] +
                      _Z[x][ y + 1][ z - 1] +
                      _Z[x][ y - 1][ z - 1];
                val = val / 8;
            }
            
            file << val << " ";
        }
        file << std::endl;
    }

    file.close();
    std::cout << "Данные EX сохранены в файле: " << filename << std::endl;
}

void Field::saveExToFileXY(const std::string &filename)
{
    std::ofstream file(filename);

    if (!file.is_open())
    {
        std::cerr << "Ошибка открытия файла для записи: " << filename << std::endl;
        return;
    }
    int x = len_x / slice;
    for (size_t y = 1; y < len_y-1; y++)
    {
        for (size_t z = 1; z < len_z-1; z++)
        {
    
    
            float val = _Y[x][y][z];
            
            if (val == 0)
            {
                val = _Y[x][y][z-1] +
                      _Y[x][y-1][ z] +
                      _Y[x][ y][ z + 1] +
                      _Y[x][ y + 1][ z] +
                      _Y[x][ y + 1][ z + 1] +
                      _Y[x][ y - 1][ z + 1] +
                      _Y[x][ y + 1][ z - 1] +
                      _Y[x][ y - 1][ z - 1];
                val = val / 8;
            }
            
            file << val << " ";
        }
        file << std::endl;
    }

    file.close();
    std::cout << "Данные EX сохранены в файле: " << filename << std::endl;
}


#endif //FLEX

#ifndef FLEX
void Field::saveExToFileZX(const std::string &filename)
{
    std::ofstream file(filename);

    if (!file.is_open())
    {
        std::cerr << "Ошибка открытия файла для записи: " << filename << std::endl;
        return;
    }
    int z = len_z / slice;
    for (size_t y = 0; y < len_y; y++)
    {
        for (size_t x = 0; x < len_x; x++)

        {
            file << _X[x][y][z] << " ";
        }
        file << std::endl;
    }
    file.close();
    std::cout << "Данные EX сохранены в файле: " << filename << std::endl;
}

void Field::saveExToFileZY(const std::string &filename)
{
    std::ofstream file(filename);

    if (!file.is_open())
    {
        std::cerr << "Ошибка открытия файла для записи: " << filename << std::endl;
        return;
    }
    int z = len_z / slice;
    for (size_t y = 0; y < len_y; y++)
    {
        for (size_t x = 0; x < len_x; x++)

        {
            file << _Y[x][y][z] << " ";
        }
        file << std::endl;
    }

    file.close();
    std::cout << "Данные EY сохранены в файле: " << filename << std::endl;
}

void Field::saveExToFileYX(const std::string &filename)
{
    std::ofstream file(filename);

    if (!file.is_open())
    {
        std::cerr << "Ошибка открытия файла для записи: " << filename << std::endl;
        return;
    }
    int y = len_y / slice;
    for (size_t x = 0; x < len_x; x++)
    {
    for (size_t z = 0; z < len_z; z++)
        {
            file << _X[x][y][z] << " ";
        }
        file << std::endl;
    }

    file.close();
    std::cout << "Данные EX сохранены в файле: " << filename << std::endl;
}

void Field::saveExToFileYZ(const std::string &filename)
{
    std::ofstream file(filename);

    if (!file.is_open())
    {
        std::cerr << "Ошибка открытия файла для записи: " << filename << std::endl;
        return;
    }
    int y = len_y / slice;
    for (size_t x = 0; x < len_x; x++)
    {
        for (size_t z = 0; z < len_z; z++)
        {
            file << _Z[x][y][z] << " ";
        }
        file << std::endl;
    }

    file.close();
    std::cout << "Данные EZ сохранены в файле: " << filename << std::endl;
}

void Field::saveExToFileXZ(const std::string &filename)
{
    std::ofstream file(filename);

    if (!file.is_open())
    {
        std::cerr << "Ошибка открытия файла для записи: " << filename << std::endl;
        return;
    }
    int x = len_x / slice;
    for (size_t y = 0; y < len_y; y++)
    {
    for (size_t z = 0; z < len_z; z++)
        {
            file << _Z[x][y][z] << " ";
        }
        file << std::endl;
    }

    file.close();
    std::cout << "Данные EX сохранены в файле: " << filename << std::endl;
}

void Field::saveExToFileXY(const std::string &filename)
{
    std::ofstream file(filename);

    if (!file.is_open())
    {
        std::cerr << "Ошибка открытия файла для записи: " << filename << std::endl;
        return;
    }
    int x = len_x / slice;
    for (size_t y = 0; y < len_y; y++)
    {
        for (size_t z = 0; z < len_z; z++)
        {
            file << _Y[x][y][z] << " ";
        }
        file << std::endl;
    }

    file.close();
    std::cout << "Данные EX сохранены в файле: " << filename << std::endl;
}

#endif
void FDTD::update(float time)
{
    time += _time;
    while (time > _time)
    {
        H.updateH(E, AH);
        E.updateE(H, J, AE, BE);
#ifdef MUR
        std::cout << "MUR\n";
        if (_time >= dt)
        {
            Mur();
        }
        GetBorderValuesMur();
#endif

#ifdef PML
    #ifdef MUR
        #error
    #else
        std::cout << "PML\n";
        SetBorderValuesPML();
        Pml();
        GetBorderValuesPML();
    #endif    
#endif
        _time += dt;
        std::cout << "time " << _time << "\n";
    }
}

void FDTD::initCoeffi()
{
for (size_t z = 0; z < len_z; z++)
    {
        for (size_t y = 0; y < len_y; y++)
    {
            for (size_t x = 0; x < len_x; x++)

            {
                AH[x][y][z] = dt / (Mu._X[x][y][z] * mu0);

                AE[x][y][z] = ((Eps._X[x][y][z] * eps0) - (0.5 * Sigm._X[x][y][z] * dt)) / ((Eps._X[x][y][z] * eps0) + (0.5 * Sigm._X[x][y][z] * dt));

                BE[x][y][z] = dt / ((Eps._X[x][y][z] * eps0) + (0.5 * Sigm._X[x][y][z] * dt));
            }
        }
    }
}
