**Data**
[[Rockfish Recruitment & Ecosystem Assessment Survey]]

**Species**
[[Pacific sanddab]] 
- Access to both YOY and adult, only using YOY from RREAS
	- No length measurements before 2011

[[Pacific hake]]
- Access to only YOY from RREAS
	- Length measurements from 1994

**Models**
- Simple GAMs
	- Either used Tweedie or Gaussian
	- Either with or without variable coefficient terms
	- PDO, NPGO, or ONI
	- Formulations:
		lncatch ~ f(year) + s(lon, lat) + s(depth) + s(DOY)
		lncatch ~ f(year) + s(lon, lat) + s(depth) + s(DOY) + s(climate index)
		lncatch ~ f(year) + s(lon, lat) + s(depth) + s(DOY) + s(lat, lon, by = climate index)
	- Results
		- [[Pacific sanddab]]
			- Gaussian: Model with s(ONI), $R^2$ = 0.356, deviance = 37.5%
			- Tweedie: Model with s(lat, lon, by = NPGO), deviance = 59.9%
		- [[Pacific hake]]
			- Gaussian: Model with s(lat, lon, by = NPGO), $R^2$ = 0.232, deviance = 25.8%
			- Tweedie: Model with s(lat, lon, by = NPGO), deviance = 47.3%
- Length binned GAMs
	- Split the data into two length bins
		- Used histograms to look at number to split at
	- Compare split models to base model using all the data
		- Had to limit due to lack of length data in some years
	- Results
		- [[Pacific sanddab]]
			- Data from 2011 onward
			- Base model: s(NPGO), deviance = 57.7%
			- Lower model: s(NPGO), deviance = 58.3%
			- Upper model: s(NPGO), deviance = 60.4%
		- [[Pacific hake]]
			- Data from 1994 onward
			- Base model: s(lat, lon, by = NPGO), deviance = 46.1%
			- Lower bin model: s(lat, lon, by = NPGO), deviance = 54.7%
			- Upper bin model: s(lat, lon, by = NPGO), deviance = 48.5%
- Notes
	- [[Pacific sanddab]]
		- Lack of inclusion of VC term
		- Improved slightly by binning lengths
	- [[Pacific hake]]
		- Improved by binning lengths


[[West Coast Groundfish Bottom Trawl Survey]]

Cleaned data to isolate species of interest
Length/age data for hake, sanddab, shortbelly, widow but not for anchovy
- Not a lot of it, quite a bit missing

**Models**
- Simple GAMs
	- Either used Tweedie or Gaussian
	- Either with or without variable coefficient terms
	- PDO, NPGO, or ONI
	- Formulations:
		lncatch ~ f(year) + s(lon, lat) + s(depth) + s(DOY) + s(temp)
		lncatch ~ f(year) + s(lon, lat) + s(depth) + s(DOY) + s(temp) + s(climate index)
		lncatch ~ f(year) + s(lon, lat) + s(depth) + s(DOY) + s(temp) + s(lat, lon, by = climate index)
	- Results
		- [[Pacific Hake]]
			- Gaussian: Model with s(lat, lon, by = NPGO), $R^2$ = 0.216, deviance = 22.5%
			- Tweedie: Model with s(lat, lon, by = NPGO), deviance = 45.5%
		- [[Pacific Sanddab]]
			- Gaussian: Model with s(lat, lon, by = NPGO), $R^2$ = 0.338, deviance = 35.1%
			- Tweedie: Model with s(lat, lon, by = NPGO), deviance = 47.3%
		- [[Anchovy]]
			- Gaussian: Model with s(lat, lon, by = NPGO), $R^2$ = 0.23, deviance = 35.9%
			- Tweedie: Model with s(lat, lon, by = NPGO), deviance = 47.3%
		- [[Shortbelly Rockfish]]
			- Gaussian: Model with s(lat, lon, by = NPGO), $R^2$ = 0.161, deviance = 19.5%
			- Tweedie: Model with s(lat, lon, by = NPGO), deviance = 47.3%
		- [[Widow Rockfish]]
			- Gaussian: Model with s(lat, lon, by = ONI), $R^2$ = 0.113, deviance = 18.5%
			- Tweedie: Model with s(lat, lon, by = NPGO), deviance = 47.3%
			
- Length binned GAMs
	- Split the data into two length bins
		- Used histograms to look at number to split at
	- Compare split models to base model using all the data
		- Had to limit due to lack of length data in some years
	- Results
		- [[Pacific hake]]
			- Base model: s(lat, lon, by = NPGO), deviance = 45.5%
			- Lower bin model: s(lat, lon, by = NPGO), deviance = 45.8%
			- Upper bin model: s(lat, lon, by = NPGO), deviance = 52.3%
		- [[Pacific sanddab]]
			- Base model: s(NPGO), deviance = 57.7%
			- Lower model: s(NPGO), deviance = 58.3%
			- Upper model: s(NPGO), deviance = 60.4%

- Notes
	- [[Pacific sanddab]]
		- Lack of inclusion of VC term
		- Improved slightly by binning lengths
	- [[Pacific hake]]
		- Improved by binning lengths