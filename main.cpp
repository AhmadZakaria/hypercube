#include <iostream>

namespace helper {
//  Windows
#ifdef _WIN32
#include <Windows.h>
double get_wall_time() {
    LARGE_INTEGER time, freq;
    if (!QueryPerformanceFrequency(&freq)) {
        //  Handle error
        return 0;
    }
    if (!QueryPerformanceCounter(&time)) {
        //  Handle error
        return 0;
    }
    return (double)time.QuadPart / freq.QuadPart;
}
double get_cpu_time() {
    FILETIME a, b, c, d;
    if (GetProcessTimes(GetCurrentProcess(), &a, &b, &c, &d) != 0) {
        //  Returns total user time.
        //  Can be tweaked to include kernel times as well.
        return
            (double)(d.dwLowDateTime |
                     ((unsigned long long)d.dwHighDateTime << 32)) * 0.0000001;
    } else {
        //  Handle error
        return 0;
    }
}

//  Posix/Linux
#else
#include <time.h>
#include <sys/time.h>
double get_wall_time() {
    struct timeval time;
    if (gettimeofday(&time, NULL)) {
        //  Handle error
        return 0;
    }
    return (double)time.tv_sec + (double)time.tv_usec * .000001;
}
double get_cpu_time() {
    return (double)clock() / CLOCKS_PER_SEC;
}
#endif
}


using namespace std;

#define P 10000
#define d 1000
#define REAL float
double start_time, end_time;

REAL get_length(REAL vec[]);
REAL get_angle(REAL a[], REAL b[], REAL lenA, REAL lenB);
void calc_distances(REAL **pts, REAL **dists_out);
REAL inline calc_expected(REAL *pts, size_t length) ;
REAL inline calc_expected(REAL **pts, size_t length) ;
void write_to_file(string filename, REAL * pts, int length);
void write_to_file(string filename, REAL ** pts, int length);
void calc_hist(REAL *pts, size_t length, string filename) ;
void calc_hist(REAL **pts, size_t length, string filename) ;


void inline start_section(string header) {
    #pragma omp master
    {
        cout << header;
        start_time = helper::get_wall_time();
    }
}

void inline end_section() {
    #pragma omp master
    {
        end_time = helper::get_wall_time();
        cout << "Done! " << end_time - start_time << " s." << endl;
    }
    #pragma omp barrier

}

int main() {

    /* initialize random seed: */
    srand (0); // for reproducability
    //srand (time(NULL));// for random

    start_section("Allocating: ");
    REAL** points = new REAL*[P];
    REAL* lengths = new REAL[P];
    REAL* angles = new REAL[P];
    REAL** dists = new REAL*[P];
    #pragma omp parallel
    {
        #pragma omp for
        for(int i = 0; i < P; ++i)
            points[i] = new REAL[d];
        #pragma omp for
        for(int i = 0; i < P; ++i)
            dists[i] = new REAL[P];

        end_section();

        start_section("Initializing: ");

        #pragma omp for
        for (int pnt = 0; pnt < P; ++pnt) {
            for (int dim = 0; dim < d; ++dim) {
                points[pnt][dim] = ((double) rand() / (RAND_MAX));
            }
        }
        end_section();

        start_section("Calculating Lengths: ");
        for (int pnt = 0; pnt < P; ++pnt) {
            lengths[pnt] = get_length(points[pnt]);
        }
        end_section();

        start_section("Calculating Angles: ");
        REAL diagonal[d];
        std::fill_n(diagonal, d, 1);
        REAL lenDiag = get_length(diagonal);

        #pragma omp for
        for (int pnt = 0; pnt < P; ++pnt) {
            angles[pnt] = get_angle(points[pnt], diagonal, lengths[pnt], lenDiag);
        }
        end_section();

        start_section("Calculating Distances: ");
        calc_distances(points, dists);
        end_section();


        start_section("Calculating Expected Values: ");
        REAL exp_length = calc_expected(lengths, P);
        REAL exp_angle = calc_expected(angles, P);
        REAL exp_dist = calc_expected(dists, P);
        end_section();

        start_section("Writing Histograms: ");

        #pragma omp single
        {
            calc_hist(lengths, P, "lengths_hist.dat");
            calc_hist(angles, P, "angles_hist.dat");
            calc_hist(dists, P, "dists_hist.dat");
        }
        end_section();

        #pragma omp single
        {
            cout << "\nExpected length: " << exp_length << endl;
            cout << "Expected angle: " << exp_angle << endl;
            cout << "Expected dist: " << exp_dist << endl <<endl;
        }
//        start_section("Writing to files: ");
//        #pragma omp single
//        {
//            write_to_file("lengths.dat", lengths, P);
//            write_to_file("angles.dat", angles, P);
//            write_to_file("dists.dat", dists, P);
//        }
//        end_section();


        start_section("Cleaning up: ");
        #pragma omp for
        for(int i = 0; i < P; ++i) {
            delete [] points[i];
            delete [] dists[i];
        }


    }
    delete [] points;
    delete [] lengths;
    delete [] angles;
    delete [] dists;

    end_section();
    return 0;
}
#include <cmath>
#include <numeric>
#include <iostream>
#include <vector>
#include <iterator>

REAL get_length(REAL vec[]) {
    return sqrt(std::inner_product(vec, vec + d, vec, 0.0));
}

REAL inline get_angle(REAL a[], REAL b[], REAL lenA, REAL lenB) {
    return acos(std::inner_product(a, a + d, b, 0.0) / (lenA * lenB));
}

REAL get_dist(REAL a[], REAL b[]) {
    REAL ret = 0.0;
#pragma unroll(d)
    for (int var = 0; var < d; ++var) {
        double dist = (*a++) - (*b++);
        ret += dist * dist;
    }
    return sqrt(ret);
}

/**
 * @brief calc_distances calculates pairwise distances by
 * calculating only the upper triangle of the distance matrix,
 * and copying it into the lower triangle
 * @param pts
 * @param dists_out
 */
void calc_distances(REAL **pts, REAL **dists_out) {
    #pragma omp for
    for (int i = 0; i < P; i++) {
        for (int j = 0; j <  i; j++) {
            REAL dis = get_dist(pts[i], pts[j]);
            dists_out[i][j] = dis;
            dists_out[j][i] = dis;
        }
    }
}

REAL inline calc_expected(REAL *pts, size_t length) {

    REAL res = 0.0;
    for (int i = 0; i < length; ++i) {
        res += pts[i];
    }
    return res / length;
}

REAL inline calc_expected(REAL **pts, size_t length) {

    REAL res = 0.0;
    for (int i = 0; i < length; ++i) {
        for (int j = 0; j < length; ++j) {
            res += pts[i][j];
        }
    }
    return res / (length * length);
}

#include <fstream>

void write_to_file(string filename, REAL * pts, int length) {
    ofstream out_file (filename);
    if (out_file.is_open()) {
        for (int i = 0; i < length; ++i) {
            out_file << pts[i] << endl;
        }
        out_file.close();
    } else cout << "Unable to open file";

}

void write_to_file(string filename, int * pts, int length) {
    ofstream out_file (filename);
    if (out_file.is_open()) {
        for (int i = 0; i < length; ++i) {
            out_file << pts[i] << endl;
        }
        out_file.close();
    } else cout << "Unable to open file";

}

void write_to_file(string filename, REAL ** pts, int length) {
    ofstream out_file (filename);
    if (out_file.is_open()) {
        for (int i = 0; i < length; ++i) {
            for (int j = 0; j < length; ++j) {
                out_file << pts[i][j] << endl;
            }
        }
        out_file.close();
    } else cout << "Unable to open file";

}
#include <float.h>

void calc_hist(REAL *pts, size_t length, string filename) {
    const int num_bins = 200;
    int counts[num_bins] = { 0 };
    float val_min = FLT_MAX, val_max = FLT_MIN;
    for (int i = 0; i != length; ++i) {
        if((pts[i]) < val_min) {
            val_min = (pts[i]);
        }
        if((pts[i]) > val_max) {
            val_max = (pts[i]);
        }
    }
    REAL bin_width = (val_max - val_min) / num_bins;
    for (int i = 0; i != length; ++i) {
        size_t bin_idx = (int)(((pts[i]) - val_min) / bin_width);
        counts[bin_idx]++;
    }
    ofstream out_file (filename);
    if (out_file.is_open()) {
        for (int i = 0; i < num_bins; ++i) {
            out_file << val_min + (i*bin_width)<<"\t"<< counts[i] << endl;
        }
        out_file.close();
    } else cout << "Unable to open file";

}

void calc_hist(REAL **pts, size_t length, string filename) {
    const int num_bins = 200;
    int counts[num_bins] = { 0 };
    float val_min = FLT_MAX, val_max = FLT_MIN;
    for (int i = 0; i < length; ++i)
        for (int j = 0; j < length; ++j) {
            if((pts[i][j]) < val_min) {
                val_min = (pts[i][j]);
            }
            if((pts[i][j]) > val_max) {
                val_max = (pts[i][j]);
            }
        }

    REAL bin_width = (val_max - val_min) / num_bins;
    for (int i = 0; i < length; ++i)
        for (int j = 0; j < length; ++j) {
            size_t bin_idx = (int)(((pts[i][j]) - val_min) / bin_width);
            counts[bin_idx]++;
        }

    ofstream out_file (filename);
    if (out_file.is_open()) {
        for (int i = 0; i < num_bins; ++i) {
            out_file << val_min + (i*bin_width)<<"\t"<< counts[i] << endl;
        }
        out_file.close();
    } else cout << "Unable to open file";}
