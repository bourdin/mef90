/* FORD MathJax configuration (see ford.md: mathjax_config).
   This file is loaded BEFORE MathJax itself, so configuration must go
   through the window.MathJax pre-load variable, not MathJax.Hub.Config.
   Enable $...$ for inline math in addition to the MathJax default \(...\). */
window.MathJax = {
    tex2jax: {
        inlineMath: [['$', '$'], ['\\(', '\\)']],
        processEscapes: true,
        /* MathJax skips pre/code by default, which leaves math in the
           colorized source listings raw; skip only the tags where
           typesetting can never make sense. */
        skipTags: ['script', 'noscript', 'style', 'textarea']
    }
};
