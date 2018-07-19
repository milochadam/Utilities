#include <algorithm>
#include <cstddef>
template <class T>
class Array2D {
    T* _array;

   public:
    const std::size_t width;
    const std::size_t height;
    Array2D(size_t x, size_t y) : width(x), height(y) {
        _array = new T[x * y];
    }

    Array2D(const Array2D& o) : width(o.width), height(o.height) {
        _array = new T[width * height];
        std::copy(o._array, o._array + o.width * o.height, _array);
    }

    ~Array2D() {
        delete[] _array;
    }

    void init() {
        for (size_t i = 0; i < width * height; ++i) {
            _array[i] = 0;
        }
    }

    T& operator()(size_t x, size_t y) {
        return _array[height * y + x];
    }

    T& operator=(const T& o) {
        if (o.width != width || o.height != height) {
            delete[] _array;
            _array = new T[o.width * o.height];
        }
        this->width = o.width;
        this->height = o.height;
        std::copy(o._array, o._array + o.width * o.height, _array);
    }
};
