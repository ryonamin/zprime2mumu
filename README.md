~ README ~
=============================================
zprime2mumu Analyzer
by R. Yonamine and F. Zenoni
=============================================

This analyzer makes use of the 'MakeClass' ROOT function.
To run the current code, do the following:

1. root -l
   Open your local root release

2. .L ZprimeLoop.C++
   Compile the macro

3. ZprimeLoop m("insert_path_of_the_ROOT_file_here")
   Create an object precising the path of the IIHE ROOT file as an argument (won't work with other TTrees)

4. m.Loop()
   Execute Loop() function.

5. Enjoy your output (if any)!
