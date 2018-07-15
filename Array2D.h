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

    ~Array2D() {
        delete[] _array;
    }

    T& operator()(size_t x, size_t y) {
        return _array[width * y + x];
    }

    void init() {
        for (size_t i = 0; i < width * height; ++i) {
            _array[i] = 0;
        }
    }
};
