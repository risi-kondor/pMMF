### pMMF

A high performance, open source, parallel MMF library written in C++.<br>
Authors:  [Risi Kondor](http://people.cs.uchicago.edu/~risi/), [Nedelina Teneva](http://people.cs.uchicago.edu/~nteneva/) and [Pramod K. Mudrakarta]().

Documentation: [soon]<br>
Gallery: [[gallery]](http://people.cs.uchicago.edu/~risi/MMF/gallery.html)<br>
BibTeX entry: [[bib]](http://people.cs.uchicago.edu/~risi/MMF/pMMF.bib)

## Dependencies

Compiling pMMF requires a [C++11](https://en.wikipedia.org/wiki/C%2B%2B11) compatible compiler such as [clang](http://clang.llvm.org). [GCC](https://gcc.gnu.org) compatibility is not guaranteed, since GCC's C++11 support is still experimental. pMMF makes extensive use of the [C++ Standard Template Library](https://en.wikipedia.org/wiki/Standard_Template_Library), which is bundled with most C++ compilers. 

Certain features of pMMF also require the following:<br>
&nbsp;&nbsp;&nbsp; [Eigen](http://eigen.tuxfamily.org) (for matrix inversion and certain tests)<br>
&nbsp;&nbsp;&nbsp; [MatIO](http://sourceforge.net/projects/matio/) (for [MATLAB](http://www.mathworks.com/products/matlab/) format input/output)<br>


## License 

pMMF is shared under the terms of the [GNU Public License V.3.0](http://www.gnu.org/licenses/gpl-3.0.en.html). The copyright to pMMF is retained by the authors, Risi Kondor, Nedelina Teneva and Pramod K. Mudrakarta, and the authors reserve the right to separately license the code, in part or in whole, to commercial entities. 
