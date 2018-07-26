//
//  timer.hpp
//  
//
//  Created by John Stanco on 3/28/18.
//
//

#ifndef timer_hpp
#define timer_hpp

#include <stdio.h>
#include <time.h>

class timer {
private:
  clock_t time;
  clock_t dt;
public:
  timer() : time(0), dt(0) {}
  void start() { dt = clock(); }
  void stop() { dt = clock() - dt; time += dt; dt = 0; }
  void reset() { time = 0; }
  void restart() { time = 0; dt = clock(); }
  double t_sec() { return time/(double)CLOCKS_PER_SEC; }
  double t_min() { return time/(double)CLOCKS_PER_SEC/60; }
  double t_hrs() { return time/(double)CLOCKS_PER_SEC/3600; }
  void print() {
    double t = time / (double)CLOCKS_PER_SEC;
    int hrs = t / 3600;
    int min = t / 60 - hrs * 60;
    double sec = t - (hrs * 3600 + min * 60);
    
    printf("hrs\t\tmin\t\tsec\n");
    printf("%d\t\t%d\t\t%.7f\n", hrs, min, sec);
  }
};

#endif /* timer_hpp */
