#ifndef HEADER_OBJ
#define HEADER_OBJ
class Obj
{
private:
    int obj_x = 75;
    int obj_y = 75;
    int obj_z = 75;
    int _obj_len_x = 50;
    int _obj_len_y = 50;
    int _obj_len_z = 50;

public:
    inline int get_obj_x() const { return obj_x; }
    inline int get_obj_y() const { return obj_y; }
    inline int get_obj_z() const { return obj_z; }
    inline int get__obj_len_x() const { return _obj_len_x; }
    inline int get__obj_len_y() const { return _obj_len_y; }
    inline int get__obj_len_z() const { return _obj_len_z; }

    void set_obj_x(const int &obj_x)
    {
        this->obj_x = obj_x;
    }
    void set_obj_y(const int &obj_y)
    {
        this->obj_y = obj_y;
    }
    void set_obj_z(const int &obj_z)
    {
        this->obj_z = obj_z;
    }
    void set_obj_len_x(const int &_obj_len_x)
    {
        this->_obj_len_x = _obj_len_x;
    }
    void set_obj_len_y(const int &_obj_len_y)
    {
        this->_obj_len_y = _obj_len_y;
    }
    void set_obj_len_z(const int &_obj_len_z)
    {
        this->_obj_len_z = _obj_len_z;
    }

public:
    Obj(int obj_x,
        int obj_y,
        int obj_z,
        int _obj_len_x,
        int _obj_len_y,
        int _obj_len_z)
    {
        this->obj_x = obj_x;
        this->obj_y = obj_y;
        this->obj_z = obj_z;
        this->_obj_len_x = _obj_len_x;
        this->_obj_len_y = _obj_len_y;
        this->_obj_len_z = _obj_len_z;
    }
};

#endif