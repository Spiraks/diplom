
#include "FDTD.h"
#include <cmath>
// #include "omp.h"

float Field::gt(size_t x, size_t y, size_t z, size_t comp)
{
    return nodes[x][y][z][comp];
}
float Field::gtO(size_t x, size_t y, size_t z, size_t comp)
{
    return nodes_1[x][y][z][comp];
}
void Field::copyNodes()
{
    nodes_1 = nodes;
};

void Field::setNode(size_t x, size_t y, size_t z, size_t comp, float value)
{
    nodes[x][y][z][comp] = value;
}

void Field::updateE(Field &H, Field &J, Field &c1, Field &c2)
{

    float k = (dt / eps0);

    // #pragma omp target teams distribute parallel for collapse(2) map(tofrom: H, J, c1,c2)
    for (size_t z = 2; z < len_z - 2; z++)
    {
        for (size_t y = 2; y < len_y - 2; y++)
        {
            for (size_t x = 2; x < len_x - 2; x++)
            {
                // Не считаем обьект
                // if (x > obj_x && x < obj_x + _obj_len_x)
                //     if (y > obj_y && y < obj_y + _obj_len_y)
                //         continue;
                // if(x>len_x - 4 && z>len_z - 4 && y>len_y - 4)
                // std::cout<<"Ex"<< x <<"-" << y << "-" << z<< "\n";
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

void Field::updateH(Field &E, Field &c0)
{

    // #pragma omp target teams distribute parallel for collapse(2) map(tofrom: E, c)
    for (size_t z = 2; z < len_z - 4; z++)
    {
        for (size_t y = 2; y < len_y - 4; y++)
        {
            for (size_t x = 2; x < len_x - 4; x++)
            {

                // Не считаем обьект
                // if (x > obj_x && x < obj_x + _obj_len_x)
                //     if (y > obj_y && y < obj_y + _obj_len_y)
                //         continue;

                this->setNode(x, y + 1, z + 1, X, this->gt(x, y + 1, z + 1, X) + c0.gt(x, y + 1, z + 1, X) * ((E.gt(x, y + 1, z + 2, Y) - E.gt(x, y + 1, z, Y)) / dz - (E.gt(x, y + 2, z + 1, Z) - E.gt(x, y, z + 1, Z)) / dy));

                this->setNode(x + 1, y, z + 1, Y, this->gt(x + 1, y, z + 1, Y) + c0.gt(x + 1, y, z + 1, Y) * ((E.gt(x + 2, y, z + 1, Z) - E.gt(x, y, z + 1, Z)) / dx - (E.gt(x + 1, y, z + 2, X) - E.gt(x + 1, y, z, X)) / dz));

                this->setNode(x + 1, y + 1, z, Z, this->gt(x + 1, y + 1, z, Z) + c0.gt(x + 1, y + 1, z, Z) * ((E.gt(x + 1, y + 2, z, X) - E.gt(x + 1, y, z, X)) / dy - (E.gt(x + 2, y + 1, z, Y) - E.gt(x, y + 1, z, Y)) / dx));
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
    for (size_t y = 1; y < len_y - 1; y++)
    {
        for (size_t x = 1; x < len_x - 1; x++)
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

void Field::saveExToFileY(const std::string &filename, int comp)
{
    std::ofstream file(filename);

    if (!file.is_open())
    {
        std::cerr << "Ошибка открытия файла для записи: " << filename << std::endl;
        return;
    }
    int y = len_y / 2;
    for (size_t x = 1; x < len_x - 1; x++)
    {
        for (size_t z = 1; z < len_z - 1; z++)
        {
            float val = gt(x, y, z, comp);
            if (val == 0)
            {
                val = gt(x - 1, y, z, comp) +
                      gt(x, y, z - 1, comp) +
                      gt(x + 1, y, z, comp) +
                      gt(x, y, z + 1, comp) +
                      gt(x + 1, y, z + 1, comp) +
                      gt(x + 1, y, z - 1, comp) +
                      gt(x - 1, y, z + 1, comp) +
                      gt(x - 1, y, z - 1, comp);
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
    for (size_t y = 0; y < len_y; y++)
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

    for (size_t y = 0; y < len_y; y++)
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
    for (size_t z = 1; z < len_z - 1; z++)
    {
        for (size_t y = 1; y < len_y - 1; y++)
        {
            float val = gt(x, y, z, comp);
            // std::cout<< y << "-" << z<< "\n";

            if (val == 0)
            {
                val = gt(x, y, z - 1, comp) +
                      gt(x, y - 1, z, comp) +
                      gt(x, y, z + 1, comp) +
                      gt(x, y + 1, z, comp) +
                      gt(x, y + 1, z + 1, comp) +
                      gt(x, y - 1, z + 1, comp) +
                      gt(x, y + 1, z - 1, comp) +
                      gt(x, y - 1, z - 1, comp);
                val = val / 8;
            }
            file << val << " ";
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
        H.updateH(E, c0);
        E.updateE(H, J, c1, c2);
#ifdef MUR
        std::cout << "MUR\n";
        if (_time >= dt)
        {
            Mur();
        }
        GetBorderValuesMur();
#endif
        // (x1 - x2) ^ 2 + (y1 - y2) ^ 2;

        _time += dt;
        std::cout << "time " << _time << "\n";
    }
}

void FDTD::initCoeffi()
{
    for (size_t z = 1; z < len_z - 2; z++)
    {
        for (size_t y = 1; y < len_y - 2; y++)
        {
            for (size_t x = 1; x < len_x - 2; x++)
            {
                for (size_t comp = X; comp <= Z; comp++)
                {
                    c0.setNode(x, y, z, comp, dt / (Mu.gt(x, y, z, comp) * mu0));
                    c1.setNode(x, y, z, comp, ((Eps.gt(x, y, z, comp) * eps0) - (0.5 * Sigm.gt(x, y, z, comp) * dt)) / ((Eps.gt(x, y, z, comp) * eps0) + (0.5 * Sigm.gt(x, y, z, comp) * dt)));
                    c2.setNode(x, y, z, comp, dt / ((Eps.gt(x, y, z, comp) * eps0) + (0.5 * Sigm.gt(x, y, z, comp) * dt)));
                }
            }
        }
    }
}

void Field::saveToFile(const std::string &filename)
{
    std::ofstream file_x(filename + "X.txt");
    std::ofstream file_y(filename + "Y.txt");
    std::ofstream file_z(filename + "Z.txt");

    if (!file_x.is_open())
    {
        std::cerr << "Ошибка открытия файла для записи: X" << filename << std::endl;
        return;
    }
    if (!file_y.is_open())
    {
        std::cerr << "Ошибка открытия файла для записи: Y" << filename << std::endl;
        return;
    }
    if (!file_z.is_open())
    {
        std::cerr << "Ошибка открытия файла для записи: Z" << filename << std::endl;
        return;
    }
    for (size_t z = 1; z < len_z - 1; z++)
    {
        for (size_t y = 1; y < len_y - 1; y++)
        {
            for (size_t x = 1; x < len_x - 1; x++)
            {

                float val = gt(x, y, z, X);
                if (val == 0)
                {
                    val = gt(x - 1, y, z, X) +
                          gt(x, y - 1, z, X) +
                          gt(x + 1, y, z, X) +
                          gt(x, y + 1, z, X) +
                          gt(x + 1, y + 1, z, X) +
                          gt(x + 1, y - 1, z, X) +
                          gt(x - 1, y + 1, z, X) +
                          gt(x - 1, y - 1, z, X);
                    val = val / 8;
                }
                file_x << val << " ";

                val = gt(x, y, z, Y);
                if (val == 0)
                {
                    val = gt(x - 1, y, z, Y) +
                          gt(x, y - 1, z, Y) +
                          gt(x + 1, y, z, Y) +
                          gt(x, y + 1, z, Y) +
                          gt(x + 1, y + 1, z, Y) +
                          gt(x + 1, y - 1, z, Y) +
                          gt(x - 1, y + 1, z, Y) +
                          gt(x - 1, y - 1, z, Y);
                    val = val / 8;
                }
                file_y << val << " ";

                val = gt(x, y, z, Z);
                if (val == 0)
                {
                    val = gt(x - 1, y, z, Z) +
                          gt(x, y - 1, z, Z) +
                          gt(x + 1, y, z, Z) +
                          gt(x, y + 1, z, Z) +
                          gt(x + 1, y + 1, z, Z) +
                          gt(x + 1, y - 1, z, Z) +
                          gt(x - 1, y + 1, z, Z) +
                          gt(x - 1, y - 1, z, Z);
                    val = val / 8;
                }
                file_z << val << " ";
            }
            file_x << std::endl;
            file_y << std::endl;
            file_z << std::endl;
        }
        file_x << std::endl;
        file_y << std::endl;
        file_z << std::endl;
    }

    file_x << std::endl;
    file_y << std::endl;
    file_z << std::endl;
    std::cout << "Данные сохранены в файле: " << filename << std::endl;
}