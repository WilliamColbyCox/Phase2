//Good Programming Practices: https://www.browserstack.com/guide/coding-standards-best-practices
//Sources:
//removing file extension from file name: https://stackoverflow.com/questions/6417817/easy-way-to-remove-extension-from-a-filename
//converting strings to ints and doubles: https://www.geeksforgeeks.org/cpp-program-for-char-to-int-conversion/
//making random unique numbers: https://www.youtube.com/watch?v=wnjYD_euSJQ
//extracting numbers from string: https://codescracker.com/cpp/program/cpp-extract-numbers-from-string.htm
//dynamic 2d arrays: https://www.geeksforgeeks.org/how-to-declare-a-2d-array-dynamically-in-c-using-new-operator/
//remove the space from nextLine: https://stackoverflow.com/questions/23834624/remove-first-and-last-character-c
//new random number algorithm: https://en.cppreference.com/w/cpp/numeric/random/uniform_int_distribution

#include <iostream>
#include <fstream>
#include <string>
#include <stdbool.h>
#include <random>
#include <vector>

class Point {
private:
    int dimensionality = 1;
    double* pointData;

    void stringToDoubleStar(std::string string) {
        double stringDouble = 0;
        pointData = new double[dimensionality];
        for (int i = 0; i < dimensionality; i++) {
            stringDouble = std::stod(string);
            pointData[i] = std::stod(string);
            //have to remove entered data point from nextLine
            if (string.find(" ") != std::string::npos) {
                string = string.substr(string.find(" "));
                string.erase(0, 1);
            }
        }
    }

    int findDimensionality(std::string string) {
        for (int i = 0; i < string.length(); i++) {
            if (string[i] == ' ') {
                dimensionality++;
            }
        }
        return dimensionality;
    }

public:
    Point(std::string pointString) {
        findDimensionality(pointString);
        stringToDoubleStar(pointString);
    }

    int getDimensionality() {
        return dimensionality;
    }

    double getValue(int dimension) {
        return pointData[dimension];
    }
};

class Cluster {
private:
    int clusterCount;
    std::vector<double> centroid;
    std::vector<Point> points;
public:
    Cluster(int clusterCount, Point newCentroid) {
        this->clusterCount = clusterCount;
        for (int i = 0; i < newCentroid.getDimensionality(); i++) {
            centroid.push_back(newCentroid.getValue(i));
        }
    }

    void addPoint(Point point) {
        points.push_back(point);
    }

    void removeAllPoints() {
        points.clear();
    }

    Point getPoint(int pointCount) {
        return points[pointCount];
    }

    double getPointValue(int pointCount, int dimension) {
        return points[pointCount].getValue(dimension);
        
    }

    int getNumberOfPointsInCluster() {
        return points.size();
    }

    void setCentroid(std::vector<double> newCentroid) {
        centroid = newCentroid;
    }

    double getCentroidValue(int dimension) {
        return centroid[dimension];
    }
};

class KMeans {
private:
    int numberOfClusters;
    int iterations;
    double convergenceThreshold;
    int runs;
    int dimensionality;
    int numberOfPoints;
    std::string outFileName;
public:
    std::vector<Cluster> clusters;
    std::vector<Point> points;

    KMeans(int numberOfClusters, int iterations, double convergenceThreshold, int runs, std::vector<Point> points, std::string outFileName) {
        this->numberOfClusters = numberOfClusters;
        this->iterations = iterations;
        this->runs = runs;
        this->points = points;
        this->convergenceThreshold = convergenceThreshold;
        this->outFileName = outFileName;
    }

    void doKMeans() {
        int iteration = 1;
        int bestRun;
        double bestSSE = 0;
        double oldSSE = 2;
        double newSSE = 1;
        numberOfPoints = points.size();
        dimensionality = points[0].getDimensionality();

        std::ofstream outfile(outFileName);
        for (int i = 1; i <= runs; i++) {
            if (i == 2) {
                bestSSE = newSSE;
            }
            outfile << "Run " + std::to_string(i) + "\n-----\n";
            initClusters();
            iteration = 1;
            while (iteration <= iterations && (oldSSE - newSSE) / oldSSE > convergenceThreshold) {
                clusterPoints();
                newSSE = resetCentroid();
                outfile << "Iteration " + std::to_string(iteration) + ": SSE = " + std::to_string(newSSE) + "\n";
                if (newSSE < bestSSE) {
                    bestSSE = newSSE;
                    bestRun = iteration;
                }
                oldSSE = newSSE;
                iteration++;
            }
        }
        outfile << "Best Run: " + std::to_string(bestRun) + ": SSE = " + std::to_string(bestSSE);
        outfile.close();
    }

    void initClusters() {
        int* selectedCenters = new int[numberOfClusters];

        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<> distrib(1, numberOfPoints);

        bool unique;
        int newRandom;
        for (int i = 0; i < numberOfClusters; i++) {
            do {
                newRandom = distrib(gen);
                unique = true;
                for (int j = 0; j < i; j++)
                    if (selectedCenters[j] == newRandom) {
                        unique = false;
                    }
            } while (!unique);
            selectedCenters[i] = newRandom;
            Cluster cluster(i, points[newRandom]);
            clusters.push_back(cluster);
        }
    }

    void clusterPoints() {
        int nearestCluster;
        for (int i = 0; i < numberOfPoints; i++) {
            nearestCluster = getNearestCluster(points[i]);
            clusters[nearestCluster].addPoint(points[i]);
        }
    }

    int getNearestCluster(Point point) {
        double distance;
        double minimumDistance;
        double sum = 0.0;
        double x;
        int nearestCluster = 0;
        for (int i = 0; i < dimensionality; i++) {
            x = clusters[0].getCentroidValue(i) - point.getValue(i);
            x = x * x;
            sum += x;
        }
        minimumDistance = sqrt(sum);
        for (int i = 1; i < numberOfClusters; i++) {
            sum = 0.0;
            for (int j = 0; j < dimensionality; j++) {
                x = clusters[i].getCentroidValue(j) - point.getValue(j);
                x = x * x;
                sum += x;
            }
            distance = sqrt(sum);
            if (minimumDistance > distance) {
                minimumDistance = distance;
                nearestCluster = i;
            }
        }
        return nearestCluster;
    }

    double resetCentroid() {
        double SSE = 0;
        for (int i = 0; i < numberOfClusters; i++) {
            int numberOfPointsInCluster = clusters[i].getNumberOfPointsInCluster();
            std::vector<double> newCentroid;
            for (int j = 0; j < dimensionality; j++) {
                double sum = 0.0;
                if (numberOfPointsInCluster > 0) {
                    for (int k = 0; k < numberOfPointsInCluster; k++) {
                        sum += clusters[i].getPoint(k).getValue(j);
                    }
                    newCentroid.push_back(sum / numberOfPointsInCluster);
                    SSE += sum / numberOfPointsInCluster;
                }
            }
            clusters[i].setCentroid(newCentroid);
        }
        return SSE;
    }
};

int main(int argc, char* argv[])
{
    //initialize arguments
    std::string fileName = argv[1];
    int clusters = std::stoi(argv[2]);
    int iterations = std::stoi(argv[3]);
    double convergenceThreshold = std::atof(argv[4]);
    int runs = std::stoi(argv[5]);

    //open file and get point info
    std::ifstream infile(fileName);
    std::string firstLine;
    getline(infile, firstLine);
    size_t lastindex = fileName.find_last_of(".");
    std::string outFileName = fileName.substr(0, lastindex) + "_output.txt";

    //get numberOfPoints and dimensionality
    int numberOfPoints = std::stoi(firstLine);
    int dimensionality = std::stoi(firstLine.substr(firstLine.find(" ")));

    //get individual data points
    std::string nextLine;
    std::vector<Point> points;
    for (int i = 0; i < numberOfPoints; i++) {
        std::getline(infile, nextLine);
        Point point(nextLine);
        points.push_back(point);
    }
    infile.close();
    KMeans kmeans(clusters, iterations, convergenceThreshold, runs, points, outFileName);
    kmeans.doKMeans();
}