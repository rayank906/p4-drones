// Project Identifier: 1761414855B69983BD8035097EFBD312EB2517F0
#include <iostream>
#include <string>
#include <vector>
#include <limits>
#include <queue>
#include <cmath>
#include <numeric>
#include <iomanip>
#include <getopt.h>
using namespace std;

enum Campus {
    Medical,
    Main
};

struct Coordinate {
    int x;
    int y;
    Campus campus;
};

struct Prims {
    Coordinate coord;
    bool visited = false;
    double dist = numeric_limits<double>::infinity();
    size_t predecessor_idx = -1;
};

class Drones {
 private:
    char mode;
    int numCoords;
    double bestDist; // optTSP
    double pathCost;
    vector<size_t> tour_idx;
    vector<size_t> path;
    vector<Coordinate> coordinates;
    vector<Prims> primCoords;
    vector<Prims> optTSPCoords;

    double mst_dist_helper(const Coordinate& a, const Coordinate& b) {
        if (a.campus != b.campus) {
            bool a_border = (a.x <= 0 && a.y <= 0 && (a.x == 0 || a.y == 0));
            bool b_border = (b.x <= 0 && b.y <= 0 && (b.x == 0 || b.y == 0));

            if (!a_border && !b_border) {
                return numeric_limits<double>::infinity();
            }
        }

        double euclid_x = a.x - b.x;
        double euclid_y = a.y - b.y;

        return (euclid_x * euclid_x) + (euclid_y * euclid_y);
    } // mst_dist_helper()

    double tsp_dist_helper(const Coordinate& a, const Coordinate& b) {
        double euclid_x = a.x - b.x;
        double euclid_y = a.y - b.y;

        return (euclid_x * euclid_x) + (euclid_y * euclid_y);
    } // tsp_dist_helper()

    double mst_opt_helper(vector<size_t>& nodes) {
        if (nodes.size() <= 1) {
            return 0;
        }
        // populate Prims
        primCoords.resize(nodes.size());
        for (size_t i = 0; i < nodes.size(); i++) {
            primCoords[i].coord = coordinates[nodes[i]];
            primCoords[i].dist = numeric_limits<double>::infinity();
            primCoords[i].visited = false;
        }

        // run Prims
        primCoords[0].dist = 0;
        for (size_t i = 0; i < primCoords.size() - 1; i++){
            size_t min_idx = -1;
            double min_dist = numeric_limits<double>::infinity();
            for (size_t j = 0; j < primCoords.size(); j++) {
                if (!primCoords[j].visited && primCoords[j].dist < min_dist) {
                    min_idx = j;
                    min_dist = primCoords[j].dist;
                }
            }
            primCoords[min_idx].visited = true;
            for (size_t k = 0; k < primCoords.size(); k++) {
                if (!primCoords[k].visited) {
                    double distance = tsp_dist_helper(primCoords[min_idx].coord, primCoords[k].coord);
                    if (distance < primCoords[k].dist) {
                        primCoords[k].dist = distance;
                        primCoords[k].predecessor_idx = min_idx;
                    }
                }
            }
        }

        double total = 0;
        for (size_t i = 1; i < primCoords.size(); i++) {
            total += sqrt(primCoords[i].dist);
        }

        return total;
    } // mst_opt_helper()

    double lower_bound_helper(size_t permLength) {
        if (pathCost >= bestDist) {
            return numeric_limits<double>::infinity();
        }
        // mst unvisited
        vector<size_t> unvisited;
        unvisited.reserve(optTSPCoords.size());
        for (size_t i = 0; i < optTSPCoords.size(); i++) {
            if (optTSPCoords[i].visited == false) {
                unvisited.push_back(i);
            }
        }

        double mstCost = 0;
        mstCost = mst_opt_helper(unvisited);

        // cost of connnecting
        double minStart = numeric_limits<double>::infinity();
        double minEnd = numeric_limits<double>::infinity();

        size_t lastNode = (permLength > 0) ? path[permLength - 1] : 0;

        for (size_t i = 0; i < unvisited.size(); i++) {
            double distFromStart = sqrt(tsp_dist_helper(coordinates[0], coordinates[unvisited[i]]));
            if (distFromStart < minStart) {
                minStart = distFromStart;
            }

            double distFromEnd = sqrt(tsp_dist_helper(coordinates[unvisited[i]], coordinates[lastNode]));
            if (distFromEnd < minEnd) {
                minEnd = distFromEnd;
            }
        }

        if (unvisited.empty()) {
            return pathCost + sqrt(tsp_dist_helper(coordinates[lastNode], coordinates[0]));
        }

        return pathCost + mstCost + minStart + minEnd;
    } // lower_bound_helper()

    bool promising(size_t permLength) {
        optTSPCoords[0].visited = true;
        for (size_t i = 0; i < permLength; i++) {
            optTSPCoords[path[i]].visited = true;
        }

        double lowerB = lower_bound_helper(permLength);

        optTSPCoords[0].visited = false;
        for (size_t i = 0; i < permLength; i++) {
            optTSPCoords[path[i]].visited = false;
        }
        return lowerB < bestDist;
    } // promising()

    void genPerms(size_t permLength) {
        if (permLength == path.size()) {
            double tourCost = pathCost + sqrt(tsp_dist_helper(coordinates[path[path.size() - 1]], coordinates[0]));

            if (tourCost < bestDist) {
                bestDist = tourCost;
                tour_idx = path;
            }
            return;
        }  // if ..complete path

        if (!promising(permLength)) {
            return;
        }  // if ..not promising

        for (size_t i = permLength; i < path.size(); ++i) {
            swap(path[permLength], path[i]);

            size_t lastNode = (permLength > 0) ? path[permLength - 1] : 0;
            size_t firstNode = path[permLength];
            double edgeCost = sqrt(tsp_dist_helper(coordinates[lastNode], coordinates[firstNode]));
            pathCost += edgeCost;

            genPerms(permLength + 1);

            pathCost -= edgeCost;

            swap(path[permLength], path[i]);
        }  // for ..unpermuted elements
    }  // genPerms()

 public:
    Drones() {
        mode = 'a';
        numCoords = 0;
        bestDist = 0.0;
        pathCost = 0.0;
    } // Drones()

    // Read in CLI
    void printHelp(char *command) {
        cout << "Usage: " << command << "-m {MST|FASTTSP|OPTTSP} | -h\n";
        cout << "P4: Drones\n";
    }  // printHelp()

    void getOptions(int argc, char **argv) {
        // These are used with getopt_long()
        opterr = static_cast<int>(false);  // Let us handle all error output for command line options
        int choice = 0;
        int index = 0;

        // NOLINTBEGIN: getopt predates C++ style, this usage is from `man getopt`
        option longOptions[] = {
            {"help", no_argument, nullptr, 'h'},
            {"mode", required_argument, nullptr, 'm'},
            {nullptr, 0, nullptr, '\0'},
        };  // longOptions[]
        // NOLINTEND

        // Fill in the double quotes, to match the mode and help options.
        while ((choice = getopt_long(argc, argv, "hm:", static_cast<option *>(longOptions), &index)) != -1) {
            switch (choice) {
            case 'h':
                printHelp(*argv);
                exit(0);

            case 'm': {
                string arg { optarg };
                if (arg == "MST") {
                    mode = 'm';
                }
                else if (arg == "FASTTSP") {
                    mode = 'f';
                }
                else if (arg == "OPTTSP") {
                    mode = 'o';
                }
                else {
                    cerr << "Error: Invalid / missing mode\n";
                    exit(1);
                }
                break;
            }

            default:
                cerr << "Error: Invalid command line option\n" << flush;
                exit(1);

            }  // switch ..choice
        }  // while
    }  // getOptions()

    void readCoords() {
        cin >> numCoords;
        coordinates.resize(numCoords);

        for (size_t i = 0; i < coordinates.size(); i++) {
            cin >> coordinates[i].x >> coordinates[i].y;
            if (coordinates[i].x < 0 && coordinates[i].y < 0) {
                coordinates[i].campus = Campus::Medical;
            }
            else {
                coordinates[i].campus = Campus::Main;
            }
        }
    } // readCoords()

    void runMode() {
        if (mode == 'm') {
            mstMode();
        }
        else if (mode == 'f') {
            fastTSPMode();
        }
        else if (mode == 'o') {
            optTSPMode();
        }
    } // runMode()

    void mstMode() {
        // populate Prims
        primCoords.resize(numCoords);
        for (size_t i = 0; i < coordinates.size(); i++) {
            primCoords[i].coord = coordinates[i];
        }

        // run Prims
        primCoords[0].dist = 0;
        for (size_t i = 0; i < primCoords.size() - 1; i++){
            size_t min_idx = -1;
            double min_dist = numeric_limits<double>::infinity();
            for (size_t j = 0; j < primCoords.size(); j++) {
                if (!primCoords[j].visited && primCoords[j].dist < min_dist) {
                    min_idx = j;
                    min_dist = primCoords[j].dist;
                }
            }
            
            // check for mst mistake
            if (min_dist == numeric_limits<double>::infinity()) {
                cerr << "Cannot construct MST\n";
                exit(1);
            }

            primCoords[min_idx].visited = true;
            for (size_t k = 0; k < primCoords.size(); k++) {
                if (!primCoords[k].visited) {
                    double distance = mst_dist_helper(primCoords[min_idx].coord, primCoords[k].coord);
                    if (distance < primCoords[k].dist) {
                        primCoords[k].dist = distance;
                        primCoords[k].predecessor_idx = min_idx;
                    }
                }
            }
        }

        // output mst
        double total = 0;
        for (size_t i = 1; i < primCoords.size(); i++) {
            if (primCoords[i].dist == numeric_limits<double>::infinity()) {
                cerr << "Cannot construct MST\n";
                exit(1);
            }
            total += sqrt(primCoords[i].dist);
        }

        cout << total << '\n';

        for (size_t i = 1; i < primCoords.size(); i++) {
            size_t idx_a = primCoords[i].predecessor_idx;
            size_t idx_b = i;

            if (idx_a > idx_b) {
                swap(idx_a, idx_b);
            }

            cout << idx_a << " " << idx_b << "\n";
        }

        primCoords = {};
    } // mstMode()

    void fastTSPMode() {
        tour_idx = {0, 1, 0};
        for (size_t k = 2; k < coordinates.size(); k++) {
            size_t best = 0;
            double min_increase = numeric_limits<double>::infinity();

            for (size_t p = 0; p < tour_idx.size() - 1; p++) {
                size_t i = tour_idx[p];
                size_t j = tour_idx[p + 1];

                double increase = sqrt(tsp_dist_helper(coordinates[i], coordinates[k])) + sqrt(tsp_dist_helper(coordinates[k], coordinates[j])) - sqrt(tsp_dist_helper(coordinates[i], coordinates[j]));

                if (increase < min_increase) {
                    min_increase = increase;
                    best = p + 1;
                }
            }
            tour_idx.insert(tour_idx.begin() + best, k);
        }
        
        // output
        double total = 0;
        for (size_t i = 0; i < tour_idx.size(); i++) {
            total += sqrt(tsp_dist_helper(coordinates[tour_idx[i]], coordinates[tour_idx[(i + 1) % tour_idx.size()]]));
        }

        bestDist = total; // for optTSP

        if (mode == 'f') {
            cout << total << '\n';

            for (size_t i = 0; i < tour_idx.size() - 1; i++) {
                cout << tour_idx[i] << " ";
            }
        }
    } // fastTSPMode()

    void optTSPMode() {
        // populate TSP
        optTSPCoords.resize(numCoords);
        for (size_t i = 0; i < coordinates.size(); i++) {
            optTSPCoords[i].coord = coordinates[i];
        }

        // init best sol
        fastTSPMode();

        // fix output to match TSP
        if (!tour_idx.empty() && tour_idx.back() == 0) {
            tour_idx.pop_back();
        }
        if (!tour_idx.empty() && tour_idx.front() == 0) {
            tour_idx.erase(tour_idx.begin());
        }

        // populate path
        path.resize(coordinates.size() - 1);
        iota(path.begin(), path.end(), 1);
        
        pathCost = 0.0;
        genPerms(0);

        // output
        cout << bestDist << '\n';
        cout << "0" << " ";
        for (size_t i = 0; i < tour_idx.size(); i++) {
            cout << tour_idx[i] << " ";
        }
        cout << "\n";
    } // optTSPMode()
};

int main(int argc, char *argv[]) {
    // This should be in all of your projects, speeds up I/O
    ios_base::sync_with_stdio(false);

    cout << std::setprecision(2); //Always show 2 decimal places
    cout << std::fixed; //Disable scientific notation for large numbers

    // init Drones
    Drones drones;

    // read in CLI
    drones.getOptions(argc, argv);

    // read coords
    drones.readCoords();

    // run mode
    drones.runMode();

    return 0;
} // main()