#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <chrono>
#include <assert.h>

using namespace std;
using namespace std::chrono;

// A cloth is defined by its width and length.
// The cloth can be rotated and checked whether it is square.
struct cloth {
    int w, h;

    // A cloth is greater than another if it is wider.
    bool operator>(const cloth& s) const {
        return (w != s.w ? w > s.w : h > s.h);
    }
};

// Cartesian coordinates of a 2-dimensional point.
// It is said that a cloth is in a x-y position if its top-left position
// is x-y.
struct point {
    int x, y;
};

// Represents the top-left point and the bottom-right point of a rectangle.
struct square {
    point p1, p2;
};

int N, W;
vector<cloth> cloths;
vector<int> widths;

// Start the chronometer.
auto start = high_resolution_clock::now();


// Read an instance and fill the global vector "cloths" with sorted values.
void read_instance(const char* file){
    ifstream in(file);
    in >> W >> N;

    cloths = vector<cloth>(N);
    widths = vector<int>(W+1, -1);

    // While reading make sure to save the input in the vector in a concrete
    // way. That is, keeping the width as the largest dimension whenever
    // possible.
    int n, w, h, it = 0;
    while(in >> n >> w >> h) {
        if (w < h and h <= W)
            swap(w, h);

        for(int i = 0; i < n; ++i)
            cloths[it++] = {w, h};
    }

    // Sort the vector according to longest width.
    sort(cloths.begin(), cloths.end(), greater<cloth>());

    // We will use the vector 'widths' as a shortcut to find the cloths we are
    // looking for faster on the vector 'cloths'.
    widths[cloths[0].w] = 0;
    for (int i = 1; i < N; ++i) {
        int curr = cloths[i].w;
        if (curr != cloths[i-1].w) widths[curr] = i;
    }
    for (int i = W; i > 0; --i)
        if (widths[i] == -1) {
            int j = i-1;
            while (j > 0 and widths[j] == -1) j--;
            widths[i] = widths[j];
        }
}

// Function that uses both global vectors to determine an efficient upper bound
// to the problem while storing the solution in a vector.
int greedy(vector<square>& solution) {
    // Local variables.
    int it, count = 0;
    int top = 0, bot = 0, right = 0;
    vector<bool> used(N, false);

    while (count < N) {
        // Place the largest cloth not used yet.
        it = 0;
        while(used[it]) it++;

        used[it] = true;
        bot += cloths[it].h;
        right = cloths[it].w;

        solution.push_back({{0, top}, {right-1, bot-1}});
        ++count;

        // Continue adding coths to the right of the previous one placed.
        // It must be fulfilled that they have less height than the first.
        it = widths[W - right];
        while (right < W and it != -1) {

            // Find the largest cloth that fits.
            while (it < N and (used[it] or cloths[it].h > bot - top)) it++;
            if (it == N) break;

            used[it] = true;
            solution.push_back({{right, top},
                                {right + cloths[it].w - 1, top + cloths[it].h - 1}});

            right += cloths[it].w;
            it = widths[W - right];
            ++count;
        }

        // Prepare next iteration.
        top = bot;
        right = 0;
    }
    // Get the longitude of the solution.
    return bot;
}


// Write the solution of the problem in the given file.
void print(int L, vector<square>& solution, const char* file){
    // Stop the chronometer and get the execution time.
    auto stop = high_resolution_clock::now();
    double time = duration_cast<milliseconds>(stop - start).count();

    ofstream out(file);

    // Convert the Wall time to seconds and round it to one decimal.
    out << static_cast<int>(time * 0.01 + 0.5) / 10.0 << endl;
    out << L << endl;
    for (auto sol : solution){
        out << sol.p1.x << " " << sol.p1.y << "   ";
        out << sol.p2.x << " " << sol.p2.y << endl;
    }
    out.close();
}


int main(int argc, char** argv){
    // Input treatment.
    if (argc == 1){
        cout << "Usage: " << argv[0] << " INPUT_FILE OUTPUT_FILE" << endl;
        exit(0);
    }
    assert(argc == 3);

    // Read input.
    read_instance(argv[1]);

    // Call greedy function.
    vector<square> solution;
    int L = greedy(solution);

    // Print the solution.
    print(L, solution, argv[2]);
}
