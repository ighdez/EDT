# EDT

EDT is a Python-based tool to construct D-efficient designs for Discrete Choice Experiments. EDT combines enough flexibility to construct from simple 2-alternative designs with few attributes, to more complex settings that may involve conditions between attributes.

While EDT designs are based on the widely-used D-efficiency criterion (see Kuhfeld, 2005), it differs from other free- or open-source efficient design tools (such as *idefix* for R) on the use of a **Random Swapping Algorithm** based on the work of Quan, Rose, Collins and Bliemer (2011), obtaining significant speed improvements to reach an optimal design, to a level that competes with well-known paid software such as NGene.

The main features of EDT are:

+ Allows to customize each attribute in terms of:
  + Attribute Levels
  + Continuous or Dummy coding (Effects coding is work-in-progress)
  + Assignement of prior parameters
  + Attribute names

+ Designs with constraints: EDT allows to define conditions over different attribute levels.
+ Designs with blocks.
+ Designs with alternative-specific constants (ASC).
+ Multiple stopping criteria (Fixed number of iterations, iterations without improvement or fixed time).
+ Allows to export the output design in an Excel file.

Any contributions to EDT are welcome via this Git, or to the email joseignaciohernandezh at gmail dot com. 

# References

+ Kuhfeld, W. F. (2005). Experimental design, efficiency, coding, and choice designs. *Marketing research methods in SAS: Experimental design, choice, conjoint, and graphical techniques*, 47-97.
+ Quan, W., Rose, J. M., Collins, A. T., & Bliemer, M. C. (2011). A comparison of algorithms for generating efficient choice experiments.