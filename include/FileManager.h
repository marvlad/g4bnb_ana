#ifndef FILEMANAGER_H
#define FILEMANAGER_H

#include <string>
#include <vector>

class FileManager {
public:
    // Constructor automatically detects if input is .root or .txt
    explicit FileManager(const std::string& path);

    const std::vector<std::string>& GetFiles() const { return files_; }

private:
    std::vector<std::string> files_;

    bool IsRootFile(const std::string& s) const;
    bool IsTextFile(const std::string& s) const;

    void LoadSingle(const std::string& filename);
    void LoadList(const std::string& listfile);
};

#endif
