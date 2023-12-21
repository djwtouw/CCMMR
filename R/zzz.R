.onUnload <- function (libpath)
{
    library.dynam.unload("CCMMR", libpath)
}


.onAttach <- function(libname, pkgname)
{
    package_citation = paste("\nTouw, D.J.W., Groenen, P.J.F., and Terada, Y.",
                             "(2022). Convex Clustering through MM: An",
                             "Efficient Algorithm to Perform Hierarchical",
                             "Clustering. arXiv preprint arXiv:2211.01877.\n")
    msg = paste("Thank you for using CCMMR!",
                "To acknowledge our work, please cite the paper:",
                package_citation, sep = "\n")
    packageStartupMessage(msg)
}
