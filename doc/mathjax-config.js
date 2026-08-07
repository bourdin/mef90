/* FORD MathJax configuration (see ford.md: mathjax_config).
   This file is loaded BEFORE MathJax itself, so configuration must go
   through the window.MathJax pre-load variable, not MathJax.Hub.Config.
   Enable $...$ for inline math in addition to the MathJax default \(...\). */
window.MathJax = {
    tex2jax: {
        inlineMath: [['$', '$'], ['\\(', '\\)']],
        processEscapes: true
    }
};
