#include "FileManager.h"
#include "Dk2nuAnalyzer.h"

int main() {

    FileManager fm("input.txt");
    auto files = fm.GetFiles();

    Dk2nuAnalyzer ana(files);

    ana.Run();
    ana.Save("analysis.root");

    return 0;
}
