(TeX-add-style-hook
 "AUSM_RHD_Notes"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-package-options
                     '(("geometry" "margin=1in") ("euscript" "mathscr")))
   (add-to-list 'LaTeX-verbatim-environments-local "lstlisting")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "lstinline")
   (add-to-list 'LaTeX-verbatim-macros-with-delims-local "lstinline")
   (TeX-run-style-hooks
    "latex2e"
    "article"
    "art10"
    "amsfonts"
    "amsmath"
    "amssymb"
    "xcolor"
    "amsthm"
    ""
    "atbegshi"
    "picture"
    "geometry"
    "mathrsfs"
    "euscript"
    "enumerate"
    "cleveref"
    "listings"
    "graphicx"
    "tikz"
    "wasysym")
   (TeX-add-symbols
    '("circled" 1))
   (LaTeX-add-labels
    "fig:mesh1")
   (LaTeX-add-amsthm-newtheorems
    "lemma"
    "sublemma"))
 :latex)

