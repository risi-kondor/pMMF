 

# pMMF

A high performance, open source, parallel MMF library written in C++.<br>
Authors:  [Risi Kondor](http://people.cs.uchicago.edu/~risi/), [Nedelina Teneva](http://people.cs.uchicago.edu/~nteneva/) and [Pramod K. Mudrakarta](https://www.cs.uchicago.edu/directory/pramod-kaushik-mudrakarta).

Documentation:&nbsp; [soon]<br>
Gallery:&nbsp; [[gallery]](http://people.cs.uchicago.edu/~risi/MMF/gallery.html)<br>
BibTeX entry:&nbsp; [[bib]](http://people.cs.uchicago.edu/~risi/MMF/pMMF.bib)<br>
Help/support email:&nbsp; pMMF.library at gmail.com<br>

## Dependencies

Compiling pMMF requires a [C++11](https://en.wikipedia.org/wiki/C%2B%2B11) compatible compiler such as [clang](http://clang.llvm.org). [GCC](https://gcc.gnu.org) compatibility is not guaranteed, since GCC's C++11 support is still experimental. pMMF makes extensive use of the [C++ Standard Template Library](https://en.wikipedia.org/wiki/Standard_Template_Library), which is bundled with most C++ compilers. 

Certain features of pMMF also require the following:<br>
&nbsp;&nbsp;&nbsp; [Eigen](http://eigen.tuxfamily.org) (for matrix inversion and certain tests)<br>
&nbsp;&nbsp;&nbsp; [MatIO](http://sourceforge.net/projects/matio/) (for [MATLAB](http://www.mathworks.com/products/matlab/) format input/output)<br>

## FAQ



<ul>
<li><b>What is MMF used for?</b>
<ul>
<li><b>Matrix compression/sketching</b> by providing a good quality low dimensional approximation to huge matrices, especially those arising in machine learning problems.</li>
<li><b>Preconditioning</b> to accelerate the solution of large linear systems.</li>
<li><b>Sparse approximation</b> by providing a multiresolution (wavelet) basis for A.</li>
<li><b>Data analysis</b> by helping to uncover the structure of A itself, for example, when A is the adjacency matrix of a complex network.</li>
</ul>

<li><b>What is Parallel MMF (pMMF)?</b>

pMMF is a fast concurrent algorithm for finding MMF factorizations. Empirially, on many classes of naturally occurring sparse matrices, pMMF runs in time close to linear in n.
</li>

<li><b>Where can I read more about it?</b>
<ul>
<li>R. Kondor, N. Teneva and V. Garg: <a href="http://people.cs.uchicago.edu/~risi/papers/KondorTenevaGargICML2014.pdf">Multiresolution Matrix Factorization</a> (ICML 2014) <a href="http://people.cs.uchicago.edu/~risi/papers/KondorTenevaGargICML2014.supplement.pdf">[supplement]</a>  <a href="http://techtalks.tv/talks/multiresolution-matrix-factorization/61063/">[video]</a></li>
<li>R. Kondor, N. Teneva and P. K. Mudrakarta: <a href="http://arxiv.org/abs/1507.04396">Parallel MMF: a Multiresolution Approach to Matrix Computation</a> (preprint, July 2015)</li>
</ul>
</li>
</ul>


## License 

pMMF is shared under the terms of the [GNU Public License V.3.0](http://www.gnu.org/licenses/gpl-3.0.en.html). The copyright to pMMF is retained by the authors, Risi Kondor, Nedelina Teneva and Pramod K. Mudrakarta, and the authors reserve the right to separately license the code, in part or in whole, to commercial entities. 
