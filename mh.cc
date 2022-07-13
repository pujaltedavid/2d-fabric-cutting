#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <chrono>
#include <assert.h>
#include <cstdlib>
#include <cmath>

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

    cloth& rotate() {
        swap(w, h);
        return *this;
    }
};

// A temperature is a number which starts at some point and decreases its value
// with a geometric law.
struct Temperature {
	double temp;
	double alpha;

	void update() {
		temp *= alpha;
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

// L is the actual minimum length needed for the cloth.
// MAX_L is the maximum possible length needed for the cloth.
int N, W, L, MAX_L;
const char* out_file;
vector<cloth> cloths;

// Start the chronometer.
auto START = high_resolution_clock::now();

// Read an instance and fill the global vector "cloths" with sorted values.
void read_instance(const char* file){
    ifstream in(file);
    in >> W >> N;

    cloths = vector<cloth>(N);
    MAX_L = 0;

    // Read each client rectangle command.
    int ni, w, h, it = 0;
    while (in >> ni >> w >> h) {
        // Place the cloths such that w >= h.
        if (w < h and h <= W) swap(w, h);

        for (int i = 0; i < ni; ++i){
            cloths[it++] = {w, h};
            MAX_L += max(w,h);
        }
    }
    // Sort by descending width.
    sort(cloths.begin(), cloths.end(), greater<cloth>());
}


// Print the Wall-clock time and the solution on "out_file".
void print(vector<square>& solution){
    auto stop = high_resolution_clock::now();
    // Time elapsed from the start of the program until now.
    double time = duration_cast<milliseconds>(stop - START).count();

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


vector<cloth> greedy() {
    int it, count = 0;
    int top = 0, bot = 0, right = 0;
    vector<bool> used(N, false);
    vector<int> widths(W+1, -1);
    vector<cloth> sol_cloth;

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

        sol_cloth.push_back(cloths[it]);
        ++count;

        // Continue adding coths to the right of the previous one placed.
        // It must be fulfilled that they have less height than the first.
        it = widths[W - right];
        while (right < W and it != -1) {

            // Find the largest cloth that fits.
            while (it < cloths.size() and (used[it] or cloths[it].h > bot - top)) it++;
            if (it == cloths.size()) break;

            used[it] = true;
            sol_cloth.push_back(cloths[it]);

            right += cloths[it].w;
            it = widths[W - right];
            ++count;
        }

        // Prepare next iteration.
        top = bot;
        right = 0;
    }
    return sol_cloth;
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
    for (int i = 0; i < MAX_L - c.h + 1; ++i) {
        for (int j = 0; j < W - c.w + 1; ++j) {
            if (not M[i][j] and fits(i, j, c.w, c.h, M))
                return {j, i};
        }
    }
    return {-1, -1};
}


// Calculate the length of the solution. If coord is true, also store the coordinates of
// the solution in sol_squares.
int calculate_length_coord(const vector<cloth>& sol, bool coord, vector<square>& sol_squares){
    vector<vector<bool>> M(MAX_L, vector<bool> (W, false));
    int length = 0;

    if (coord) sol_squares.clear();
    for (cloth c : sol){
        point p = place(c, M);
        assert(p.x != -1);
        length = max(length, p.y+c.h);

        if (coord) sol_squares.push_back({p, {p.x + c.w - 1, p.y + c.h - 1}});
    }
    return length;
}


// Return the probability computed by Boltzmann distribution.
double boltzmann(const int l1, const int l2, const Temperature& t){
	return exp(-(l2*1.0-l1*1.0)/t.temp);
}


void swap_solution(vector<cloth>& sol, int k, bool k_r, int m, bool m_r){
	cloth aux = sol[k];
    // Make sure that if the cloth has to be rotated, the new width
    // doesn't exceed W.
	sol[k] = (k_r and sol[m].h < W) ? sol[m].rotate() : sol[m];
	sol[m] = (m_r and aux.h < W) ? aux.rotate() : aux;
}


// Return a neighbor of the solution computed by swapping two cloths randomly.
// The cloths swapped can be rotated randomly.
vector<cloth> get_neighbor(const vector<cloth>& sol){
	int k = rand() % N;
	int m;

	do {
		m = rand() % N;
	} while (sol[k] == sol[m]);

	bool k_r = rand() % 2;
	bool m_r = rand() % 2;

	vector<cloth> newsol = sol;

	swap_solution(newsol, k, k_r, m, m_r);

	return newsol;
}


// Auxiliar function to help debugging or visualizing the solutions found.
void print_console(const vector<cloth>& newsol){
    cout << "\n\n";
	for (int i = 0; i < 75; ++i) cout << "- ";
	cout << "\nL = " << L << "\n\n";
    int i = 0;
	for (cloth c : newsol){
        cout << c.w << (c.w >= 10 ? " " : "  ") << c.h << (c.h >= 10 ? " " : "  ") << "   ";
        ++i;
        if (i%17 == 0) cout << endl;
    }
    cout << "\n\n\n";
}


// Find a solution based on Simulated Annealing metaheuristic.
void mh(vector<cloth>& sol){
    vector<square> sol_coord; // Contains the coordinates of the solution.

	print_console(sol);

    Temperature t = {3000.0, 0.95};

	while (true){
		vector<cloth> newsol = get_neighbor(sol);

		int l1 = calculate_length_coord(sol, false, sol_coord);
		int l2 = calculate_length_coord(newsol, false, sol_coord);

		if (l2 <= l1){
			// Accept the new solution.
			sol = newsol;

			if (l2 < L){
				L = l2;
                // Calculate coordinates to print in the file
                calculate_length_coord(newsol, true, sol_coord);

				print(sol_coord);
				print_console(newsol);
			}
		}
		else {
			// Accept worsening move with boltzmann probability.
			double p = boltzmann(l1, l2, t);
			double eval = rand() % 10000; // Use 4 decimals for the probability.
			if (eval < p*10000) sol = newsol;
		}
        t.update();
	}
}

int main(int argc, char** argv){
    if (argc == 1){
        cout << "Usage: " << argv[0] << " INPUT_FILE WIDTH" << endl;
        exit(0);
    }

	assert(argc == 3);

    // Set seed value for the random calls
    srand((uint) system_clock::now().time_since_epoch().count());

    read_instance(argv[1]);
    out_file = argv[2]; // Store the output filename in globals.

    // Calculate a greedy solution.
    vector<cloth> sol_greedy = greedy();
    vector<square> sol_coord;

    // Print the coordinates of the greedy solution.
    L = calculate_length_coord(sol_greedy, true, sol_coord);
	print(sol_coord);

	// Start the methaeuristic from the greedy solution.
    mh(sol_greedy);
}
