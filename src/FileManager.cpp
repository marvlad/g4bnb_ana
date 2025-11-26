#include "FileManager.h"
#include <fstream>
#include <iostream>
#include <algorithm>
#include <cctype>

FileManager::FileManager(const std::string& path) {
    if (IsRootFile(path)) {
        LoadSingle(path);
    }
    else if (IsTextFile(path)) {
        LoadList(path);
    }
    else {
        std::cerr << "ERROR: Unknown file type " << path << std::endl;
    }
}

bool FileManager::IsRootFile(const std::string& s) const {
    return s.size() > 5 && s.substr(s.size() - 5) == ".root";
}

bool FileManager::IsTextFile(const std::string& s) const {
    return s.size() > 4 && s.substr(s.size() - 4) == ".txt";
}

void FileManager::LoadSingle(const std::string& filename) {
    files_.push_back(filename);
}

void FileManager::LoadList(const std::string& listfile) {
    std::ifstream infile(listfile);
    if (!infile.is_open()) {
        std::cerr << "ERROR: Cannot open list file " << listfile << std::endl;
        return;
    }

    std::string line;
    while (std::getline(infile, line)) {
        // trim whitespace
        line.erase(std::remove_if(line.begin(), line.end(),
                [](unsigned char c){ return std::isspace(c); }),
                line.end());

        if (!line.empty())
            files_.push_back(line);
    }
}
