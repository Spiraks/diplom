
#include "FDTD.h"
#include <cmath>

float Field::getNode(float x, float y, float z, int comp)
{
    return nodes[(int)(x)][(int)(y)][(int)(z)][comp];
}
float Field::getNodeOld(float x, float y, float z, int comp)
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
                              c1.getNode(x + 1, y, z, X) * this->getNode(x + 1, y, z, X) +
                                  c2.getNode(x + 1, y, z, X) *
                                      ((H.getNode(x + 1, y + 1, z, Z) - H.getNode(x + 1, y - 1, z, Z)) / dy -
                                       (H.getNode(x + 1, y, z + 1, Y) - H.getNode(x + 1, y, z - 1, Y)) / dz) -
                                  J.getNode(x + 1, y, z, X));

                this->setNode(x, y + 1, z, Y, c1.getNode(x, y + 1, z, Y) * this->getNode(x, y + 1, z, Y) + c2.getNode(x, y + 1, z, Y) * ((H.getNode(x, y + 1, z + 1, X) - H.getNode(x, y + 1, z - 1, X)) / dz - (H.getNode(x + 1, y + 1, z, Z) - H.getNode(x - 1, y + 1, z, Z)) / dx) - J.getNode(x, y + 1, z, Y));

                this->setNode(x, y, z + 1, Z, c1.getNode(x, y, z + 1, Z) * this->getNode(x, y, z + 1, Z) + c2.getNode(x, y, z + 1, Z) * ((H.getNode(x + 1, y, z + 1, Y) - H.getNode(x - 1, y, z + 1, Y)) / dx - (H.getNode(x, y + 1, z + 1, X) - H.getNode(x, y - 1, z + 1, X)) / dy) - J.getNode(x, y, z + 1, Z));
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

                this->setNode(x, y + 1, z + 1, X, this->getNode(x, y + 1, z + 1, X) + c.getNode(x, y + 1, z + 1, X) * ((E.getNode(x, y + 1, z + 2, Y) - E.getNode(x, y + 1, z, Y)) / dz - (E.getNode(x, y + 2, z + 1, Z) - E.getNode(x, y, z + 1, Z)) / dy));

                this->setNode(x + 1, y, z + 1, Y, this->getNode(x + 1, y, z + 1, Y) + c.getNode(x + 1, y, z + 1, Y) * ((E.getNode(x + 2, y, z + 1, Z) - E.getNode(x, y, z + 1, Z)) / dx - (E.getNode(x + 1, y, z + 2, X) - E.getNode(x + 1, y, z, X)) / dz));

                this->setNode(x + 1, y + 1, z, Z, this->getNode(x + 1, y + 1, z, Z) + c.getNode(x + 1, y + 1, z, Z) * ((E.getNode(x + 1, y + 2, z, X) - E.getNode(x + 1, y, z, X)) / dy - (E.getNode(x + 2, y + 1, z, Y) - E.getNode(x, y + 1, z, Y)) / dx));
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
            float val = getNode(x, y, z, comp);
            if (val == 0)
            {
                val = getNode(x - 1, y, z, comp) +
                      getNode(x, y - 1, z, comp) +
                      getNode(x + 1, y, z, comp) +
                      getNode(x, y + 1, z, comp) +
                      getNode(x + 1, y + 1, z, comp) +
                      getNode(x + 1, y - 1, z, comp) +
                      getNode(x - 1, y + 1, z, comp) +
                      getNode(x - 1, y - 1, z, comp);
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
        file << getNode(x, y, z, comp) << " ";
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
            file << getNode(x, y, z, comp) << " ";
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
            file << getNode(x, y, z, comp) << " ";
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
                    c.setNode(x, y, z, comp, dt / (Mu.getNode(x, y, z, comp) * mu0));
                    c1.setNode(x, y, z, comp, ((Eps.getNode(x, y, z, comp) * epsilon0) - (0.5 * Sigm.getNode(x, y, z, comp) * dt)) / ((Eps.getNode(x, y, z, comp) * epsilon0) + (0.5 * Sigm.getNode(x, y, z, comp) * dt)));
                    c2.setNode(x, y, z, comp, dt / ((Eps.getNode(x, y, z, comp) * epsilon0) + (0.5 * Sigm.getNode(x, y, z, comp) * dt)));
                }
            }
        }
    }
}