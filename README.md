#  data reduction for large atmospheric satellite datasets
This repository contains R code for reducting large atmospheric satellite datasets.

Authors: Xiaoling Liu, August Weinbren, and Scot Miller. (Note: we want this code repository to be a collaborative effort, and the more authors/contributors, the better!)

----------------------------------------------------------
Introduction and overview
----------------------------------------------------------

This code repository contains R scripts for running the data reduction process. 
The overall strategy is to first characterize the decorrelation length among the satellite observation and second, use local kriging to interpolate the satellite observations to a number of locations that is smaller than the original dataset. The choice of locations is informed by the decorrelation length: we retain fewer locations in regions where the observations are correlated over longer distances and more locations in regions with a shorter decorrelation length.

Most aspects of these scripts do not need to be customized or edited by the user. However, a few aspects of the scripts need to be customized for your work, and those sections of the scripts are basically in the INITIALIZATION section and CONSTANT VARIABLES section. 

----------------------------------------------------------
Fair use
----------------------------------------------------------

Please send us an email (xiaolingliu96 [at] gmail.com, august.weinbren96 [at] gmail.com, and scot.m.miller [at] gmail.com) if you are thinking of using these scripts in your research. We would love to hear more about what you are working on and if we can help in any way. 

Also, please make sure to cite the companion article if you do use these scripts in your work.


----------------------------------------------------------
List of scripts in the repository
----------------------------------------------------------

SCRIPT: reduce.R <br>
PURPOSE: Data Reduction for large atmospheric satellite datasets. <br>
CALLED BY: None. <br>
CALLS: CalculateRange.R <br>


----------------------------------------------------------
References
----------------------------------------------------------


----------------------------------------------------------
Contact information
----------------------------------------------------------

Xiaoling Liu <br>
Johns Hopkins University <br>
Department of Environmental Health and Engineering <br>
Email: xiaolingliu96 [at] gmail.com <br>
http://greenhousegaslab.org/ 

August Weinbren <br>
Johns Hopkins University <br>
Department of Environmental Health and Engineering <br>
Email: august.weinbren96@gmail.com <br>
http://greenhousegaslab.org/ 

Scot Miller <br>
Johns Hopkins University <br>
Department of Environmental Health and Engineering <br>
Email: smill191 [at] jhu.edu OR scot.m.miller [at] gmail.com <br>
http://greenhousegaslab.org/ 


