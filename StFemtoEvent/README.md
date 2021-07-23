# StFemtoEvent

The StFemtoEvent set of classes (library) allows to store and
work with filename.femtoDst.root files.

*****************************************************************

Maintainer:   Grigory Nigmatkulov

Institution:  National Research Nuclear University MEPhI

Date:         November 15, 2018

E-mail:       nigmatkulov@gmail.com ; ganigmatkulov@mephi.ru

Contributors:
  G. Nigmatkulov, N. Ermakov, E. Alpatov, E. Khyzhniak

*****************************************************************

Brief description of how to compile and run the analysis
over femtoDst on your laptop (or computer farm):

## Compilation:
a) System has to have ROOT preinstalled (should work with versions 5 and 6).

b) There is a Makefile stored in the StFemtoEvent directory.In order to compile the codes one needs to run:

```
make
```
c) Since it is the most commonly used, the g++ compiler is used for the compilation. However, one can also switch it to clang++ (has been tested).

d) After the compilation is finished the shared library libStFemtoDst.so will be created.

e) If you have some errors during the compilation please contact the femtoDst maintainer.

## Processing femtoDst:

There are two ROOT macroses with example of how to make a simple analysis on femtoDst. They are stored in the StFemtoEvent/macros directory and called: FemtoDstAnalyzer.C and RunAnalyzer.C.

Lets assume that one has and input file(s) InputFile with a name.femtoDst.root or a list of femtoDst files, called name.lis(t), and StFemtoEvent compiled (i.e. there is libStFemtoDst.so). There are 2 possible processing scenarios depending on the ROOT version:

### ROOT 5:

One should run ROOT from the terminal:
```
[myterm]> root RunAnalyzer.C\(\"InputFile\"\)
```

Or run this macros from the interactive session:
```
[myterm]> root
root [0].x RunAnalyzer.C("InputFile")
```

Any of aforementioned ways load libStFemtoDst.so, compile and run RunAnalyzer.C. After the processing the remove dictionary and library created by ACLiC.

### ROOT 6:

Since ROOT 6 does not have CINT there is some extra flexibility on how to analyze the data. The first one is listed above, i.e. a). The second option is to run the FemtoDstAnalyzer.C macro directly.

Either from the terminal:
```
root FemtoDstAnalyzer.C\(\"InputFile"\)
```

Or from the root session:
```
[myterm]> root
root [0].x FemtoDstAnalyzer.C("InputFile")
```

For any questions or with any suggestions please contact the package maintainer.

## License
This project is licensed under the 3-Clause BSD License - see the [LICENSE.md](LICENSE.md) file for details
