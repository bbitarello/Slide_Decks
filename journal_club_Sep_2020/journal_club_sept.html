<!DOCTYPE html>
<html lang="" xml:lang="">
  <head>
    <title>Theoretical and empirical quantification of the accuracy of polygenic scores in ancestry divergent populations</title>
    <meta charset="utf-8" />
    <meta name="author" content="Bárbara Domingues Bitarello" />
    <script src="libs/header-attrs/header-attrs.js"></script>
    <link rel="stylesheet" href="xaringan-themer.css" type="text/css" />
  </head>
  <body>
    <textarea id="source">
class: center, middle, inverse, title-slide

# Theoretical and empirical quantification of the accuracy of polygenic scores in ancestry divergent populations
## Journal Club
### Bárbara Domingues Bitarello
### Perelman School of Medicine, University of Pennsylvania

---





class: left, top

# Background

* GWAS are heavily biased towards European ancestries
* Limits use of GWAS-derived PRS in other ancestries and admixed populations
--

* Duncan et al.: average PRS accuracy across traits 64% lower in individuals of African ancestry

* Martin et al. : 78% lower for African ancestry
--

* Decrease is monotonic with genetic distance (and LD differences) from EUR

* Some attempts to model this LOA (loss of accuracy) in animal breeds exist
--

* Major factors: LD, allele frequencies
* Other factors: population specific causal variants, GXE
* Quantification of relative contribution of different factors still lacking
---
class: left, top

# Overview

### Goals:

* `Theoretical understanding` of predictive capacity of trans-ancestry PRS
* Relative contributions of different factors
* Theoretical relative accuracy (RA) of trans-ancestry PRS as a function of `popgen parameters`

--

### Summary of work:

* They developed three implementations of a `deterministic formula` for RA depending on available data
* Requires `GWAS summary statistics and reference panels`
* Theory evaluation through **simulations** show RA predictor is unbiased
* Applied to **8 traits with diverse architectures**
---

# Methods
## Samples &amp; Quality Controls

--

* all the data in this study is UKBB
--

* UKBB individuals projected onto 1000G PCA (EUR, EAS, SAS, AFR) and assigned population it is closest to (3 PCs)
--

* `10,000 British used for validation` and `~10,000 Irish used for testing`. Remaining British used for GWAS
* HapMap3 SNPs present in at least 5 copies in each population ~990K
* 1KG used for reference panel for imputation in non-Europeans
* Relatedness filter: GRM, GTCA
* Final: ~ 9.4K SAS, 7K AFR, 2.2K EAS
---
class: left, top

#Methods
##Simulations

* Based on true UKBB genotypes
* Phenotypes: **Additive model** `\(y=g+e\)`
* `\(e\sim N(0,1-h^2\)`, phenotypic variance across populations is 1
* `\(M_{c}={1000,5000,10000}\)` and `\(h^2={0.25, 0.5}\)`
* `\(M_c\)` variants `sampled from HapMap3` SNPs
* `\(\rho_b=1\)`: correlation of effect sizes
* Each `\(\beta\)` sampled from `\(\sim N(0,\frac{h^2}{2p(1-p)M_c})\)`
* For each individual:

`$$g=\sum_{j=1}^{M}x_{jl}\beta_{j}$$` 

`$$x_{jl}=\frac{x_j-2p_j}{\sqrt{2p_j(1-p_j)}}$$`


---
class: left, top

#Methods
##Simulations

* GWAS and LD clumping with PLINK: `\(p&lt;5e-8\)`, `\(r^2&gt; 0.01\)`, 0.2 Mb
* Causal variants left out to mimic imperfect LD
* Pop 1: discovery; Pop 2: target/test
* Negative selection: `\(\beta\sim N(0, (2p(1-p)^2\delta^2)\)`. 
* `\(\delta^2\)` is variance of causal effect sizes
  * `\(S_1=S_2=-0.5\)`
  * `\(S_1=-0.75; S_2=-0.5\)`
  * `\(S_1=-0.5; S_2=-0.75\)`
  
* GWAS sizes: from 100k to 3000K, by 40k
* 1KG and UKBB imputed as imputation panel

---
class: left, top
#Methods: Deterministic PRS RA
* `\(RA=R^2_{2}/R^2_{1}\)`
* `Eq.1`: assumes causal variants are known **(pred1)**

`$$RA_{pred1}\approx \left[ \frac{\rho_{b}^{2} h^{2}_{2}}{h^{2}_{1}} \right]X$$`
--
`$$RA_{pred1}\approx \left[ \frac{\rho_{b}^{2} h^{2}_{2}}{h^{2}_{1}} \right] X \left( \frac{\sum_{k=1}^{M_{T}}\sqrt{\frac{p_{k,2}(1-p_{k,2})}{p_{k,1}(1-p_{k,1})}}[\sum_{j=1}^{M_{C}}r_{jk,1}r_{jk,1}]}{\sum_{k=1}^{M_{T}}(\sum_{j=1}^{M_{C}}r^{2}_{jk,1})}\right)X$$`

--
`$$RA_{pred1}\approx \left[ \frac{\rho_{b}^{2} h^{2}_{2}}{h^{2}_{1}} \right] X \left( \frac{\sum_{k=1}^{M_{T}}\sqrt{\frac{p_{k,2}(1-p_{k,2})}{p_{k,1}(1-p_{k,1})}}[\sum_{j=1}^{M_{C}}r_{jk,1}r_{jk,2}]}{\sum_{k=1}^{M_{T}}(\sum_{j=1}^{M_{C}}r^{2}_{jk,1})}\right) X \left[ \frac{var(PGS1)}{var(PGS2)}\right]$$`

---
class: left, top
#Methods: Deterministic PRS RA
* `Eq. 2`: heuristic method for candidate causal variants **(pred2)**
  * 100Kb windows, `\(r2&gt;0.45\)` with GWS SNPs
* Once candidate causal variants are identified for each GWS SNP, replace `\(r^{2}_{jk,1}\)` and `\(r_{jk,1}r_{jk,2}\)` with the average over all causal variants

`$$RA_{pred2}\approx \left[ \frac{\rho_{b}^{2} h^{2}_{2}}{h^{2}_{1}} \right] X \left( \frac{\sum_{k=1}^{M_{T}}\overline{r_{k,1}r_{k,2}}\sqrt{\frac{p_{k,2}(1-p_{k,2})}{p_{k,1}(1-p_{k,1})}}}{\sum_{k=1}^{M_{T}}\overline{r^{2}_{k,1}}}\right) X \left[ \frac{\sum_{k=1}^{M_{T}}p_{k,1}(1-p_{k,1})\hat{\beta}^{2}_{k}}{\sum_{k=1}^{M_{T}}p_{k,2}(1-p_{k,2})\hat{\beta}^{2}_{k}}\right]$$`

---
class: left, top
#Methods: Deterministic PRS RA
* `Naive`: assumes GWS SNPs are the causal SNPs **(pred3)**

`$$RA_{pred3}\approx \left[ \frac{\rho_{b}^{2} h^{2}_{2}}{h^{2}_{1}} \right] X \frac{var(PGS1)}{var(PGS2)}$$`
--

`$$RA_{pred3}\approx \left[ \frac{\rho_{b}^{2} h^{2}_{2}}{h^{2}_{1}} \right] X \frac{1}{M_{T}}\left(\sum_{k=1}^{M_{T}}\sqrt{\frac{p_{k,2}(1-p_{k,2})}{p_{k,1}(1-p_{k,1})}}\right)^2     X \frac{var(PGS1)}{var(PGS2)}$$`

* A function of the heterozygosity ratio
---
class: inverse, left, top
#Methods: Empirical PRS RA


![](table1.png)

*phenotypes regressed by sex, age, recuritment center, genotype batch, 10Pcs. 
* residuals reverse normal transformed
* outliers removed

---
class: inverse, center, middle
#Methods: Summary

![](S1.png)
---
class: inverse, center, middle
#Methods: Summary

![](S2.png)

---
class: left, top
##RESULTS
###Expected RA across ancestries
* assumptions: causal variants shared without loss of generality

--

* effect sizes allowed to vary

--

* expected `\(R^2\)` function of LD, MAF, `\(h^2\)`, `\(\rho_b\)`

--

* `\(h^2\)`, `\(\rho_b\)` require geno and phenotype from given dataset, not always available

--

* assumes  `\(h^2_{1}=h^2_{2}\)`, `\(\rho_b=1\)` **PROBABLY NOT TRUE in reality**

--

* `\(LOA=1-RA\)` 
* RA only predicts `\(LOA\)` explained by MAF and LD

--


---
class: inverse, left, top
###Performance in simulations: Accuracy

![](S3.png)

---
class: inverse, left, top
###Performance in simulations: Relative accuracy
![](Fig1.png)
---
class: inverse, left, top
###Performance in simulations: `\(h^2\)`
![](S4.png)

---
class: inverse, left, top
###Performance in simulations: RA with UKBB as reference panel
![](S5.png)

---
class: inverse, left, top
###Performance in simulations: Allelic Frequencies correlations 
![](S6.png)

---
class: inverse, left, top
###Performancein simulations: LD score correlations
* `\(LD=1+\bar{r}^{2}+m\)`; `\(m\)` is number of SNPs in 10Kb window
![](S7.png)

---
class: inverse, left, top
###Performance in simulations: Size of GWAS
![](S8.png)

---
class: inverse, left, top
##Performance in simulations: Negative selection
![](Fig2.png)
* A: `\(S_1 = S_2 =−0.5\)`
* B: `\(S_1 =−0.5, S_2 =−0.75\)`
* C: `\(S_1 =−0.75, S_2 =−0.5\)`

---
class: inverse, left, top
##Performance in real data: `Eq.2`

![](table2A.png)

---
class: inverse, left, top
##Performance in real data: `Eq.2`

![](table2B.png)


---
class: inverse, middle, top
##Performance in real data: `Eq.2`
![](fig3.png)



---
class: inverse, middle, top
##Performance in real data: more SNPs does not improve RA
![](s10.png)
---
class: left, top
## Limitations
* LD and MAF estimated from **reference panels**
* Only common SNPs, and **rare variants** more likely to be population-specific
* PGS from GWS SNPs (not whole genome)
* RA predictions might be inflated if there's strong **epistatic interactions between causal variants in LD** `or` if causal effect size is function of LD differences between populations
* Fine-mapping would be better than heuristic method in real data
--

## Take home
* RA predicitons here are `upper bounds` for how much LOA can occur due to LD and MAF alone
* large GWAS for individual ancestries would be ideal
* trans-ancestry GWAS more realistic goal: good for fine-mapping and discovering new causal/associated variants
* admixed populations: trans-ancestry GWAS with parent populations would be ideal
---

class: center, middle

# That's all folks!



    </textarea>
<style data-target="print-only">@media screen {.remark-slide-container{display:block;}.remark-slide-scaler{box-shadow:none;}}</style>
<script src="https://remarkjs.com/downloads/remark-latest.min.js"></script>
<script>var slideshow = remark.create({
"highlightStyle": "github",
"highlightLines": true,
"countIncrementalSlides": false
});
if (window.HTMLWidgets) slideshow.on('afterShowSlide', function (slide) {
  window.dispatchEvent(new Event('resize'));
});
(function(d) {
  var s = d.createElement("style"), r = d.querySelector(".remark-slide-scaler");
  if (!r) return;
  s.type = "text/css"; s.innerHTML = "@page {size: " + r.style.width + " " + r.style.height +"; }";
  d.head.appendChild(s);
})(document);

(function(d) {
  var el = d.getElementsByClassName("remark-slides-area");
  if (!el) return;
  var slide, slides = slideshow.getSlides(), els = el[0].children;
  for (var i = 1; i < slides.length; i++) {
    slide = slides[i];
    if (slide.properties.continued === "true" || slide.properties.count === "false") {
      els[i - 1].className += ' has-continuation';
    }
  }
  var s = d.createElement("style");
  s.type = "text/css"; s.innerHTML = "@media print { .has-continuation { display: none; } }";
  d.head.appendChild(s);
})(document);
// delete the temporary CSS (for displaying all slides initially) when the user
// starts to view slides
(function() {
  var deleted = false;
  slideshow.on('beforeShowSlide', function(slide) {
    if (deleted) return;
    var sheets = document.styleSheets, node;
    for (var i = 0; i < sheets.length; i++) {
      node = sheets[i].ownerNode;
      if (node.dataset["target"] !== "print-only") continue;
      node.parentNode.removeChild(node);
    }
    deleted = true;
  });
})();
(function() {
  "use strict"
  // Replace <script> tags in slides area to make them executable
  var scripts = document.querySelectorAll(
    '.remark-slides-area .remark-slide-container script'
  );
  if (!scripts.length) return;
  for (var i = 0; i < scripts.length; i++) {
    var s = document.createElement('script');
    var code = document.createTextNode(scripts[i].textContent);
    s.appendChild(code);
    var scriptAttrs = scripts[i].attributes;
    for (var j = 0; j < scriptAttrs.length; j++) {
      s.setAttribute(scriptAttrs[j].name, scriptAttrs[j].value);
    }
    scripts[i].parentElement.replaceChild(s, scripts[i]);
  }
})();
(function() {
  var links = document.getElementsByTagName('a');
  for (var i = 0; i < links.length; i++) {
    if (/^(https?:)?\/\//.test(links[i].getAttribute('href'))) {
      links[i].target = '_blank';
    }
  }
})();
// adds .remark-code-has-line-highlighted class to <pre> parent elements
// of code chunks containing highlighted lines with class .remark-code-line-highlighted
(function(d) {
  const hlines = d.querySelectorAll('.remark-code-line-highlighted');
  const preParents = [];
  const findPreParent = function(line, p = 0) {
    if (p > 1) return null; // traverse up no further than grandparent
    const el = line.parentElement;
    return el.tagName === "PRE" ? el : findPreParent(el, ++p);
  };

  for (let line of hlines) {
    let pre = findPreParent(line);
    if (pre && !preParents.includes(pre)) preParents.push(pre);
  }
  preParents.forEach(p => p.classList.add("remark-code-has-line-highlighted"));
})(document);</script>

<script>
slideshow._releaseMath = function(el) {
  var i, text, code, codes = el.getElementsByTagName('code');
  for (i = 0; i < codes.length;) {
    code = codes[i];
    if (code.parentNode.tagName !== 'PRE' && code.childElementCount === 0) {
      text = code.textContent;
      if (/^\\\((.|\s)+\\\)$/.test(text) || /^\\\[(.|\s)+\\\]$/.test(text) ||
          /^\$\$(.|\s)+\$\$$/.test(text) ||
          /^\\begin\{([^}]+)\}(.|\s)+\\end\{[^}]+\}$/.test(text)) {
        code.outerHTML = code.innerHTML;  // remove <code></code>
        continue;
      }
    }
    i++;
  }
};
slideshow._releaseMath(document);
</script>
<!-- dynamically load mathjax for compatibility with self-contained -->
<script>
(function () {
  var script = document.createElement('script');
  script.type = 'text/javascript';
  script.src  = 'https://mathjax.rstudio.com/latest/MathJax.js?config=TeX-MML-AM_CHTML';
  if (location.protocol !== 'file:' && /^https?:/.test(script.src))
    script.src  = script.src.replace(/^https?:/, '');
  document.getElementsByTagName('head')[0].appendChild(script);
})();
</script>
  </body>
</html>
