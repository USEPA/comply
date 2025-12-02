# Code Organization

The COMPLY code in this repository is organized around programming languages. There is a Fortran version, R version, and additional code that builds a user interface using R Shiny. 

# **Background & Application**

### **Background**
In 1985, the EPA asked the National Council on Radiation Protection and Measurements (NCRP) to develop simple screening methods for assessing compliance with the Clean Air Act by users of small quantities of radionuclides. NCRP published these procedures in 1986 and 1989 in Addendum to NCRP Commentary No. 3 on Determining Compliance with Environmental Standards (NCRP 89). EPA's COMPLY model was developed based on the procedures in Addendum to NCRP Commentary No. 3 on Determining Compliance with Environmental Standards.

### **Application**
COMPLY calculates the effective dose equivalent (ede) from radionuclides released from stacks and vents. Atmospheric concentrations are estimated using a Gaussian plume model and equations that account for building wake effects.

### **Strengths & Limitations**
The COMPLY computer software includes four increasing levels of complexity. A user can demonstrate compliance at any level. Level 1 requests the least amount of information: emission concentration or quantity possessed. However, under Level 1, "worst case" assumptions are used in the dose estimate. Level 4 requests the most information and uses site specific data instead of assuming the worst. The values estimated by this screening model are strictly for comparison with environmental standards and are not intended to represent actual doses to real people.

Use of COMPLY is not appropriate for estimating doses from sources that do not originate from an exhaust stack or vent.

***

Go to the [US EPA COMPLY website](https://www.epa.gov/radiation/comply) for more information.

