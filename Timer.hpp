#include <chrono>

using namespace std::chrono;
class Timer {
   public:
    Timer() {}

   private:
    high_resolution_clock::time_point start;
    high_resolution_clock::time_point stop;

   public:
    void tic() { start = high_resolution_clock::now(); }
    void toc() { stop = high_resolution_clock::now(); }
    milliseconds getTimeInMiliseconds() {
        return duration_cast<milliseconds>(stop - start);
    }
};
