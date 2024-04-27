
#include "FDTD.h"
#include <cmath>

float Field::gt(float x, float y, float z, int comp)
{
    return nodes[(int)(x)][(int)(y)][(int)(z)][comp];
}
float Field::gtO(float x, float y, float z, int comp)
{
    return nodes_old[(int)(x)][(int)(y)][(int)(z)][comp];
}
void Field::copyNodes()
{
    nodes_old = nodes;
};

void Field::setNode(float x, float y, float z, int comp, float value)
{
    nodes[(int)(x)][(int)(y)][(int)(z)][comp] = value;
}

void Field::updateE(Field &H, Field &J, Field &c1, Field &c2)
{

    float k = (dt / epsilon0);
    for (float z = 2; z < len_z - 2; z++)
    {
        for (float y = 2; y < len_y - 2; y++)
        {
            for (float x = 2; x < len_x - 2; x++)
            {
                // Не считаем обьект
                // if (x > obj_x && x < obj_x + _obj_len_x)
                //     if (y > obj_y && y < obj_y + _obj_len_y)
                //         continue;
                this->setNode(x + 1, y, z, X,
                              c1.gt(x + 1, y, z, X) * this->gt(x + 1, y, z, X) +
                                  c2.gt(x + 1, y, z, X) *
                                      ((H.gt(x + 1, y + 1, z, Z) - H.gt(x + 1, y - 1, z, Z)) / dy -
                                       (H.gt(x + 1, y, z + 1, Y) - H.gt(x + 1, y, z - 1, Y)) / dz) -
                                  J.gt(x + 1, y, z, X));

                this->setNode(x, y + 1, z, Y, c1.gt(x, y + 1, z, Y) * this->gt(x, y + 1, z, Y) + c2.gt(x, y + 1, z, Y) * ((H.gt(x, y + 1, z + 1, X) - H.gt(x, y + 1, z - 1, X)) / dz - (H.gt(x + 1, y + 1, z, Z) - H.gt(x - 1, y + 1, z, Z)) / dx) - J.gt(x, y + 1, z, Y));

                this->setNode(x, y, z + 1, Z, c1.gt(x, y, z + 1, Z) * this->gt(x, y, z + 1, Z) + c2.gt(x, y, z + 1, Z) * ((H.gt(x + 1, y, z + 1, Y) - H.gt(x - 1, y, z + 1, Y)) / dx - (H.gt(x, y + 1, z + 1, X) - H.gt(x, y - 1, z + 1, X)) / dy) - J.gt(x, y, z + 1, Z));
            }
        }
    }
}

void Field::updateH(Field &E, Field &c)
{
    for (float z = 1; z < len_z - 2; z++)
    {
        for (float y = 1; y < len_y - 2; y++)
        {
            for (float x = 1; x < len_x - 2; x++)
            {

                // Не считаем обьект
                // if (x > obj_x && x < obj_x + _obj_len_x)
                //     if (y > obj_y && y < obj_y + _obj_len_y)
                //         continue;

                this->setNode(x, y + 1, z + 1, X, this->gt(x, y + 1, z + 1, X) + c.gt(x, y + 1, z + 1, X) * ((E.gt(x, y + 1, z + 2, Y) - E.gt(x, y + 1, z, Y)) / dz - (E.gt(x, y + 2, z + 1, Z) - E.gt(x, y, z + 1, Z)) / dy));

                this->setNode(x + 1, y, z + 1, Y, this->gt(x + 1, y, z + 1, Y) + c.gt(x + 1, y, z + 1, Y) * ((E.gt(x + 2, y, z + 1, Z) - E.gt(x, y, z + 1, Z)) / dx - (E.gt(x + 1, y, z + 2, X) - E.gt(x + 1, y, z, X)) / dz));

                this->setNode(x + 1, y + 1, z, Z, this->gt(x + 1, y + 1, z, Z) + c.gt(x + 1, y + 1, z, Z) * ((E.gt(x + 1, y + 2, z, X) - E.gt(x + 1, y, z, X)) / dy - (E.gt(x + 2, y + 1, z, Y) - E.gt(x, y + 1, z, Y)) / dx));
            }
        }
    }
}

void Field::saveExToFileZ(const std::string &filename, int comp)
{
    std::ofstream file(filename);

    if (!file.is_open())
    {
        std::cerr << "Ошибка открытия файла для записи: " << filename << std::endl;
        return;
    }
    int z = len_z / 2;
    for (int y = 1; y < len_y - 1; y++)
    {
        for (int x = 1; x < len_x - 1; x++)
        {
            float val = gt(x, y, z, comp);
            if (val == 0)
            {
                val = gt(x - 1, y, z, comp) +
                      gt(x, y - 1, z, comp) +
                      gt(x + 1, y, z, comp) +
                      gt(x, y + 1, z, comp) +
                      gt(x + 1, y + 1, z, comp) +
                      gt(x + 1, y - 1, z, comp) +
                      gt(x - 1, y + 1, z, comp) +
                      gt(x - 1, y - 1, z, comp);
                val = val / 8;
            }
            file << val << " ";
        }
        file << std::endl;
    }

    file.close();
    std::cout << "Данные E" << comp << " сохранены в файле: " << filename << std::endl;
}

void Field::saveGraphZ(const std::string &filename, int comp)
{
    std::ofstream file(filename);

    if (!file.is_open())
    {
        std::cerr << "Ошибка открытия файла для записи: " << filename << std::endl;
        return;
    }
    int z = len_z / 2;
    int x = len_x / 2;
    for (int y = 0; y < len_y; y++)
    {
        file << gt(x, y, z, comp) << " ";
    }

    file.close();
    std::cout << "Данные E" << comp << " сохранены в файле: " << filename << std::endl;
}

void Field::saveRangeY()
{
    std::ofstream file("/home/alex/src/diplom_FDTD/range.txt");

    if (!file.is_open())
    {
        std::cerr << "Ошибка открытия файла для записи: "
                  << "/home/alex/src/diplom_FDTD/range.txt" << std::endl;
        return;
    }

    for (int y = 0; y < len_y; y++)
    {
        file << y << " ";
    }

    file.close();
}

void Field::saveExToFileX(const std::string &filename, int comp)
{
    std::ofstream file(filename);

    if (!file.is_open())
    {
        std::cerr << "Ошибка открытия файла для записи: " << filename << std::endl;
        return;
    }
    int x = len_x / 2;
    for (int z = 0; z < len_z; z++)
    {
        for (int y = 0; y < len_y; y++)
        {
            file << gt(x, y, z, comp) << " ";
        }
        file << std::endl;
    }

    file.close();
    std::cout << "Данные E" << comp << " сохранены в файле: " << filename << std::endl;
}

void Field::saveExToFileY(const std::string &filename, int comp)
{
    std::ofstream file(filename);

    if (!file.is_open())
    {
        std::cerr << "Ошибка открытия файла для записи: " << filename << std::endl;
        return;
    }
    int y = len_y / 2;
    for (int x = 1; x < len_x - 1; x++)
    {
        for (int z = 1; z < len_z - 1; z++)
        {
            file << gt(x, y, z, comp) << " ";
        }
        file << std::endl;
    }

    file.close();
    std::cout << "Данные E" << comp << " сохранены в файле: " << filename << std::endl;
}

void FDTD::update(float time)
{
    time += _time;
    while (time > _time)
    {
        H.updateH(E, c);
        E.updateE(H, J, c1, c2);
#ifdef MUR
        std::cout << "MUR\n";
        if (_time >= dt)
        {
            MursFirstABC();
        }
        GetBorderValues();
#endif
        // (x1 - x2) ^ 2 + (y1 - y2) ^ 2;

        _time += dt;
        std::cout << "time " << _time << "\n";
    }
}
void FDTD::initCoeffi()
{
    for (float z = 1; z < len_z - 2; z++)
    {
        for (float y = 1; y < len_y - 2; y++)
        {
            for (float x = 1; x < len_x - 2; x++)
            {
                for (size_t comp = X; comp <= Z; comp++)
                {
                    c.setNode(x, y, z, comp, dt / (Mu.gt(x, y, z, comp) * mu0));
                    c1.setNode(x, y, z, comp, ((Eps.gt(x, y, z, comp) * epsilon0) - (0.5 * Sigm.gt(x, y, z, comp) * dt)) / ((Eps.gt(x, y, z, comp) * epsilon0) + (0.5 * Sigm.gt(x, y, z, comp) * dt)));
                    c2.setNode(x, y, z, comp, dt / ((Eps.gt(x, y, z, comp) * epsilon0) + (0.5 * Sigm.gt(x, y, z, comp) * dt)));
                }
            }
        }
    }
}