#include <iostream>
#include <vector>
#include <fstream>
using namespace std;

struct point {
    int x, y;
};

struct square {
    point p1;
    point p2;
};

int L, W, N;
bool B;
vector<square> solution;

void read_instance(const char* file){
    ifstream in(file);
    double useless;
    in >> useless >> L;
    
    solution = vector<square> (N);  

    int x1, y1, x2, y2;
    int i = 0;
    while(in >> x1 >> y1 >> x2 >> y2) {
        solution[i] = {{x1,y1},{x2,y2}};
        ++i;
    }
}

void unefficient_but_useful_print(){
    // create nest
    vector< vector<char>> pict(L, vector<char>(W, '#'));

    // fill nest
    int size = solution.size();
    int wasted = W*L;
    for (int c = 0; c < size; ++c){
        int maxi = solution[c].p2.y;
        for (int i = solution[c].p1.y; i <= maxi; ++i){
            int maxj = solution[c].p2.x;
            for (int j = solution[c].p1.x; j <= maxj; ++j){
                if (B) pict[i][j] = '$'+c%86;
                else pict[i][j] = '0'+c%10;
                --wasted;
            }
        }
    }
    // print nest
    cout << "\nX: |" << (L >= 10 ? " " : "");
    for (int i = 0; i < W; ++i) cout << (i > 0 ? " " : "") << i;
    cout << " |\n___|";
    for (int i = 0; i < W*2+1; ++i) cout << "_";
    cout << "|___\n";
    for (int i = 0; i < L; ++i){
        cout << i << " " << (i<10 and L>10 ? " " : "") << "|";
        for (int j = 0; j < W; ++j){
            cout << " " << pict[i][j];
        }
        cout << " |";
        if (i == L/2){
            for(int k = 0; k < 40; ++k) cout << " ";
            cout << "|  USEFUL PRINTFUL POWER HERE TO HELP  |";
        }
        else if (i == L/2 - 1 or i == L/2 + 1){
            for(int k = 0; k < 40; ++k) cout << " ";
            for(int k = 0; k < 40; ++k) cout << "-";
            
        }
        cout << endl;
    } cout << endl;
    cout << W*L << " total cubes" << endl;
    cout << W*L - wasted << " used cubes" << endl;
    cout << wasted << " wasted cubes" << endl;
    cout << "current length        = " << L << endl;
    cout << "min impossible length = " << (W*L-wasted) / W + (W*L-wasted % W == 0 ? 0 : 1)-1 << endl;
    cout << (L == (W*L-wasted) / W + (W*L-wasted % W == 0 ? 0 : 1)-1 ? "Woah efficient" : "") << endl;
    
    cout << endl;
}

int main(int args, char** argv){
    if (args != 5){
        cout << "Usage: " << argv[0] << " INPUT_FILE WIDTH N (bool ? ascii_char : numbers" << endl;
        exit(0);
    }

    W = atoi(argv[2]);
    N = atoi(argv[3]);
    B = atoi(argv[4]);

    read_instance(argv[1]);
    unefficient_but_useful_print();
}