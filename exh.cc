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

    bool operator==(const cloth& s) const {
        return (w == s.w and h == s.h);
    }

    bool is_square() const{
        return w == h;
    }

    void rotate() {
        swap(w, h);
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

int N, W, L; // L is the actual minimum length needed for the cloth.
const char* out_file;
vector<cloth> cloths;

// Start the chronometer.
auto start = high_resolution_clock::now();

// Read an instance and fill the global vector "cloths" with sorted values.
void read_instance(const char* file){
    ifstream in(file);
    in >> W >> N;

    cloths = vector<cloth>(N);

    // Read each client rectangle command.
    int ni, w, h, it = 0;
    while (in >> ni >> w >> h) {
        // Place the cloths such that w >= h.
        if (w < h and h <= W) swap(w, h);

        for (int i = 0; i < ni; ++i)
            cloths[it++] = {w, h};
    }
    // Sort by descending width.
    sort(cloths.begin(), cloths.end(), greater<cloth>());
}

// Print the Wall-clock time and the solution on "out_file".
void print(vector<square>& solution){
    auto stop = high_resolution_clock::now();
    // Time elapsed from the start of the program until now.
    double time = duration_cast<milliseconds>(stop - start).count();

    ofstream out(out_file);

    // Convert the Wall time to seconds and round it to one decimal.
    out << static_cast<int>(time * 0.01 + 0.5) / 10.0 << endl;
    out << L << endl;

    for (auto sol : solution){
        out << sol.p1.x << " " << sol.p1.y << "   ";
        out << sol.p2.x << " " << sol.p2.y << endl;
    }
    out.close();
}

void greedy() {
    int it, count = 0;
    int top = 0, bot = 0, right = 0;
    vector<bool> used(N, false);
    vector<int> widths(W+1, -1);
    vector<square> solution;

    // The vector 'widths' is used as a shortcut to find the cloths we are
    // looking for with a specific width. It has "cloths" positions.
    widths[cloths[0].w] = 0;
    for (int i = 1; i < cloths.size(); ++i) {
        int curr = cloths[i].w;
        if (curr != cloths[i-1].w) widths[curr] = i;
    }
    // Fill the invalid positions with the previous value, so when looking
    // for a width w, if there is none, look for w-1 and so on.
    for (int i = W; i > 0; --i)
        if (widths[i] == -1) {
            int j = i-1;
            while (j > 0 and widths[j] == -1) j--;
            widths[i] = widths[j];
        }

    while (count < N) {
        // Place the largest cloth not used yet.
        it = 0;
        while (used[it]) it++;

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
            while (it < cloths.size() and (used[it] or cloths[it].h > bot - top)) it++;
            if (it == cloths.size()) break;

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
    // Save and print the greedy solution.
    L = bot;
    print(solution);
}

// "insert" = true,  insert in M the cloth in position p, defined by its
//                   width and length.
// "insert" = false, remove the cloth from M.
void replace(point p, int w, int h, vector<vector<bool> >& M, bool insert){
    for (int i = p.y; i < p.y + h; ++i)
        for (int j = p.x; j < p.x + w; ++j)
            M[i][j] = insert;
}

// Return true if the cloth fits in the i-j position. If true, insert it in M.
bool fits(int i, int j, int w, int h, vector<vector<bool>>& M) {
    // Check that the positions are not used by another cloths.
    for (int k = i; k < i + h; ++k)
        for (int z = j; z < j + w; ++z)
            if (M[k][z]) return false;

    replace({j,i}, w, h, M, true);
    return true;
}

// Place the cloth in M and return the position that has been placed,
// return a null value if can't be inserted.
point place(cloth c, vector<vector<bool>>& M) {
    for (int i = 0; i < L - c.h + 1; ++i) {
        for (int j = 0; j < W - c.w + 1; ++j) {
            if (not M[i][j] and fits(i, j, c.w, c.h, M))
                return {j, i};
        }
    }
    return {-1, -1};
}

// Return true if the i-th cloth has width "w" and height "h" or viceversa.
bool called(int i, int w, int h){
     return (w == cloths[i].w and h == cloths[i].h) or
            (w == cloths[i].h and h == cloths[i].w);
}

// Find optimal solutions using an exhaustive search.
void exhaustive (int k, int l, vector<vector<bool>>& M, vector<bool>& used,
                 vector<square>& solution) {
    if (l > L) return; // Upper bound.

    if (k == N){
        // Update L and print the current solution if it's better.
        if (l < L) {
            L = l;
            print(solution);
        }
        return;
    }
    // "w" and "h" indicates the width and height of the previous cloth
    // called recursivelly.
    int w = 0, h = 0;
    for (int i = 0; i < N; ++i){
        // If there are two or more cloths with the same dimensions, calling
        // recursively for only one of them is necessary. Since they are
        // ordered, the repeated ones are consecutive, that way it is only
        // needed to check the previous one.
        if (not used[i] and not called(i, w, h)) {
            used[i] = true;
            cloth c = cloths[i]; // To avoid multiple memory access.

            // Place the actual cloth and call the function recursively.
            point p = place(c, M);
            if (p.x != -1) {
                solution[k] = {p, {p.x + c.w - 1, p.y + c.h - 1}};
                exhaustive(k+1, max(l, p.y + c.h), M, used, solution);
                replace(p, c.w, c.h, M, false); // Remove the actual cloth.
            }

            // Do the same with the cloth rotated.
            if (not c.is_square()){
                c.rotate();
                p = place(c, M);
                if (p.x != -1) {
                    solution[k] = {p, {p.x + c.w - 1, p.y + c.h - 1}};
                    exhaustive(k+1, max(l, p.y + c.h), M, used, solution);
                    replace(p, c.w, c.h, M, false);
                }
            }

            used[i] = false;

            // Save the dimensions of the calculated cloth.
            w = c.h;
            h = c.w;
        }
    }
}


int main(int argc, char** argv){
    if (argc == 1){
        cout << "Usage: " << argv[0] << " INPUT_FILE WIDTH" << endl;
        exit(0);
    }

    assert(argc == 3);

    read_instance(argv[1]);
    out_file = argv[2]; // Store the output filename in globals.

    // Calculate a greedy solution to set a non-infinite L value to start from.
    greedy();

    vector<square> solution(N);
    // Represents an W*L roll of cloth indicating wether a position is
    // ocuped or not.
    vector<vector<bool>> M(L, vector<bool> (W, false));
    vector<bool> used(N, false);
    exhaustive(0, 0, M, used, solution);
}
