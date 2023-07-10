#' Unemployment Duration in Germany
#' @name UnempDuration
#' @docType data
#' @keywords dataset
#' @description  a cross-section of 21685 individuals from 1996 to 1997 
#' @format a tibble containing:
#' - duration: the duration of the unemployment spell in days
#' - censored: a factor with levels \code{yes} if the spell is censored, \code{no} otherwise
#' - gender: a factor with levels \code{male} and \code{female}
#' - age: the age
#' - wage: the last daily wage before unemployment
#' @source The Royal Statistical Society Datasets Website \url{http://www.blackwellpublishing.com/rss/}
#' @references
#' \insertRef{WICH:WILK:08}{micsr}
#' @importFrom Rdpack reprompt
NULL

#' Unemployment duration
#' @name Unemployment
#' @docType data
#' @keywords dataset
#' @description  a cross-section of 452 individuals from 1993 
#' @format a tibble containing:
#' - duration: duration of first spell of unemployment, t, in weeks
#' - spell: 1 if spell is complete
#' - race: one of nonwhite, white
#' - sex: one of male, female
#' - reason: reason for unemployment, one of new (new entrant), lose (job loser), leave (job leaver), reentr (labor force reentrant)
#' - search: 'yes' if (1) the unemployment spell is completed between the first and second surveys and number of methods used to search > average number of methods used across all records in the sample, or, (2) for individuals who remain unemployed for consecutive surveys, if the number of methods used is strictly nondecreasing at all survey points, and is strictly increasing at least at one survey point
#' - pubemp: 'yes' if an individual used a public employment agency to search for work at any survey points relating to the individuals first unemployment spell
#' - ftp1: 1  if an individual is searching for full time work at survey 1
#' - ftp2: 1  if an individual is searching for full time work at survey 2
#' - ftp3: 1  if an individual is searching for full time work at survey 3
#' - ftp4: 1  if an individual is searching for full time work at survey 4
#' - nobs: number of observations on the first spell of unemployment for the record
#' @source Journal of Applied Econometrics Data Archive : \url{http://qed.econ.queensu.ca/jae/}
#' @references
#' \insertRef{ROME:99}{micsr}
#' @importFrom Rdpack reprompt
NULL

#' Adoptees education 
#' @name adoptees
#' @docType data
#' @keywords dataset
#' @description  a cross-section of 16481 individuals from 1992 
#' @format a tibble containing:
#' - adopt: is the child adopted ?
#' - educ: years of education
#' - coll: college graduate
#' - stinsch: still in school
#' - gender: one of 'male' and 'female'
#' - age: the age in years
#' - feduc: father's education
#' - meduc: mother's education
#' - fcoll: is the father college graduate ?
#' - mcoll: is the mother college graduate ?
#' - linc: log of income
#' @source American Economic Association Data Archive : \url{https://www.aeaweb.org/aer/}
#' @references
#' \insertRef{PLUG:04}{micsr}
#' @importFrom Rdpack reprompt
NULL

#' Racial Discrimination in the Sharing Economy
#' @name airbnb
#' @docType data
#' @keywords dataset
#' @description  a cross-section of 6235 applications from 2015 
#' @format a tibble containing:
#' - acceptance: \code{yes} if the application is acccepted
#' - guest_race: guest race (\code{black} or \code{other}
#' - guest_gender: guest gender
#' - host_race: host race (\code{black} or \code{other}
#' - host_gender: host gender
#' - multlistings: guest has several listings
#' - shared: a dummy for shared properties
#' - tenreviews: host has ten or more reviews
#' - price: the price for a night
#' - city: the city, one of \code{Balitmore}, \code{Dallas}, \code{Los-Angeles}, \code{St-Louis} or \code{Washington}
#' @source American Economic Association Data Archive : \url{https://www.aeaweb.org/aer/}
#' @references
#' \insertRef{EDEL:LUCA:SVIR:17}{micsr}
#' @importFrom Rdpack reprompt
NULL

#' Apple production
#' @name apples
#' @docType data
#' @keywords dataset
#' @description  yearly observations of 173 farms from 1984 to 1986 
#' @format a tibble containing:
#' - id: farm's id
#' - year: year
#' - capital: capital stock
#' - labor: quantity of labor
#' - materials: quantity of materials
#' - apples: production of apples
#' - otherprod: other productions
#' - pc: price of capital
#' - pl: price of labor
#' - pm: price of materials
#' @source Journal of Applied Econometrics Data Archive : \url{http://qed.econ.queensu.ca/jae/}
#' @references
#' \insertRef{IVAL:LADO:OSSA:SIMI:96}{micsr}
#' @importFrom Rdpack reprompt
NULL

#' Economic effects of terrorism in the Basque country
#' @name basque_country
#' @docType data
#' @keywords dataset
#' @description  a pseudo-panel of 17 regions from 1955 to 1997 
#' @format a tibble containing:
#' - id: region's id
#' - region: regions's name
#' - year: the year
#' - gdpcap: gdp per capita
#' - agriculture: NA
#' - energy: 
#' - industry: NA
#' - construction: NA
#' - services: NA
#' - administration: NA
#' - illit_educ: NA
#' - prim_educ: NA
#' - medium_educ: NA
#' - high_educ: NA
#' - popdens: population density
#' @source Synth package
#' @references
#' \insertRef{ABAD:GARD:03}{micsr}
#' @importFrom Rdpack reprompt
NULL

#' Cigarette smoking and birth weight
#' @name birthwt
#' @docType data
#' @keywords dataset
#' @description  a cross-section of 1388 individuals from 1988 
#' @format a tibble containing:
#' - birthwt: birth weight
#' - cigarettes: number of cigarettes smoked per day during pregnancy
#' - parity: birth order
#' - race: a factor with levels other and white
#' - sex: a factor with levels female and male
#' - edmother: number of years of education of the mother
#' - edfather: number of years of education of the father
#' - faminc: family income
#' - cigtax: per-pack state excise tax on cigarettes
#' @source kindly provided by John Mullahy
#' @references
#' \insertRef{MULL:97}{micsr}
#' @importFrom Rdpack reprompt
NULL

#' Leniency and Cartels Enforcement
#' @name cartels
#' @docType data
#' @keywords dataset
#' @description  a time series of 40 semesterly observations from 1985 to 2005 
#' @format a tibble containing:
#' - semiyear: the semester
#' - ncaught: the number of cartels discoveries
#' - len: one if the period postdates august 10, 1993, new leniency program introduced by the Department Of Justice
#' - time: a time trend, from 1 to 40
#' - lent: a time trend starting one period after the introduction of the leniency program, zero before
#' - dgdp: semestrial growth of the gdp
#' - funds: average antitrust division budget allocation
#' - fines: total corporate fines issued by the antitrust division during the previous fiscal year
#' @source American Economic Association Data Archive : \url{https://www.aeaweb.org/aer/}
#' @references
#' \insertRef{MILL:09}{micsr}
#' @importFrom Rdpack reprompt
NULL

#' Intergenerational transmission of charitable giving
#' @name charitable
#' @docType data
#' @keywords dataset
#' @description  a cross-section of 2384 households from 2001 
#' @format a tibble containing:
#' - donation: the amount of charitable giving
#' - donparents: the amount of charitable giving of the parents
#' - education: the level of education of household's head, a factor with levels `less_high_school`, `high_school`, `some_college`, `college`, `post_college`
#' - religion: a factor with levels `none`, `catholic`, `protestant`, `jewish` and `other`.
#' - income: income
#' - married: a dummy for married couples
#' - south: a dummy for households living in the south
#' @source this data set was kindly provided by Mark Ottoni Wilhelm
#' @references
#' \insertRef{WILH:08}{micsr}
#' @importFrom Rdpack reprompt
NULL

#' Cigarette smoking behaviour
#' @name cigmales
#' @docType data
#' @keywords dataset
#' @description  a cross-section of 6160 individuals from 1979 to 1980 
#' @format a tibble containing:
#' - cigarettes: number of daily cigarettes smoked
#' - habit: smoking habit stock measure
#' - price: state-level average per-pack price of cigarettes in 1979
#' - restaurant: an indicator of whether the individual's state of residence had restrictions on smoking in restaurants in place in 1979
#' - income: family income in thousands
#' - age: age in years
#' - educ: schooling in years
#' - famsize: number of family members
#' - race: a factor with levels "other" and "white"
#' - reslgth: number of years the state's restaurant smoking restrictions had been in place in 1979
#' - lagprice: one-year lag of cigarette price
#' @source kindly provided by John Mullahy
#' @references
#' \insertRef{MULL:97}{micsr}
#' @importFrom Rdpack reprompt
NULL

#' Populist Persuasion, Evidence from Father Coughlin
#' @name coughlin
#' @docType data
#' @keywords dataset
#' @description  a cross-section of 3102 counties from 1936 
#' @format a tibble containing:
#' - state: state name
#' - county: county name
#' - id: county id
#' - fdr: share of votes for Franklin D. Roosevelt in the 1936 presidential election
#' - signal: signal strength of the Father Coughlin's radio program in 1936,
#' - signal_free: hypothetical signal strenght (on an earth free of topographic obstacles)
#' - male: share of males
#' - native_whites: share of native whites
#' - foreign_whites: share of foreign whites
#' - blacks: share of blacks
#' - urban: share of urbans
#' - elderlies: share of the population aged 65 or more
#' - catholics: share of catholics
#' - illiterate: share of illeterates
#' - unemp: unemployment rate
#' - income: occupational income score
#' - radio: share of radio owners
#' - manuf: share of manufacturing workers
#' - agri: share of agricultural workers
#' - farm_size: average farm size
#' - land_value: land value per acre
#' - tenant_acre: share of tenant acres
#' - area: area
#' - elev: mean elevation
#' - rug: rug
#' - pop: population
#' - past_democrat: average vote share of the democratic party
#' - past_republican: average vote share of the republican party
#' - past_turnout: aerage voter turnout during 1920-1928
#' @source American Economic Association Data Archive : \url{https://www.aeaweb.org/aer/}
#' @references
#' \insertRef{WANG:21}{micsr}
#' @importFrom Rdpack reprompt
NULL

#' Doctor visits in Australia
#' @name doctor_aus
#' @docType data
#' @keywords dataset
#' @description  a cross-section of 5190 individuals from 1977 to 1978 
#' @format a tibble containing:
#' - sex: a factor with levels `male` and `female`
#' - age: age in years divided by 100
#' - income: annual income in tens of thousands of dollars
#' - insurance: insurance contract, a factor with levels `medlevy`, `levyplus` (private insurance), `freepor` (free government insurance due to low-income) and `freerepa` (free government insurance due to old-age)
#' - illness: number of illness in past 2 weeks
#' - actdays: number of days of reduced activity in past 2 weeks due to illness or injury
#' - hscore: general health score using Goldberg's method (from 0 to 12)
#' - chcond: chronic condition (`np` : no problem, `nla` : chronic condition but not limited in activity, `la` : chronic condition and milited in activity)
#' - doctorco: number of consultations with a doctor or specialist in the past 2 weeks
#' - nondocco: number of consultations with non-doctor health professionals (chemist, optician, physiotherapist, social worker, district community nurse, chiropodist or chiropractor) in the past 2 weeks
#' - hospadmi: number of admissions to a hospital, psychiatric hospital, nursing or convalescent home in the past 12 months (up to 5 or more admissions which is coded as 5)
#' - hospdays: number of nights in a hospital, etc.  during most recent admission: taken, where appropriate, as the mid-point of the intervals 1, 2, 3, 4, 5, 6, 7, 8-14, 15-30, 31-60, 61-79 with 80 or more admissions coded as 80. If no admission in past 12 months then equals zero.
#' - medicine: total number of prescribed and nonprescribed medications used in past 2 days
#' - prescrib: total number of prescribed medications used in past 2 days
#' - nonpresc: total number of nonprescribed medications used in past 2 days
#' @source http://cameron.econ.ucdavis.edu/racd/racddata.html
#' @references
#' \insertRef{CAME:TRIV:86}{micsr}
#' 
#' \insertRef{CAME:JOHA:97}{micsr}
#' 
#' \insertRef{CAME:TRIV:MILN:PIGG:88}{micsr}
#' 
#' \insertRef{CAME:TRIV:13}{micsr}
#' @importFrom Rdpack reprompt
NULL

#' Physician advice on alcohol consumption
#' @name drinks
#' @docType data
#' @keywords dataset
#' @description  a cross-section of 2467 individuals from 1990 
#' @format a tibble containing:
#' - drinks: number of drinks in the past 2 weeks
#' - advice: 1 if reveived a drining advice
#' - age: age in 10 years cathegories
#' - race: a factor with levels `white`, `black` and `other`
#' - marital: marital status, one of `single`, `married`, `widow`, `separated`
#' - region: one of `west`, `northeast`, `midwest` and `south`
#' - empstatus: one of `other`, `emp` and `unemp`
#' - limits: limits on daily activities, one of `none`, `some` and `major`
#' - income: monthly income ($1000)
#' - educ: education in years
#' - medicare: insurance through medicare
#' - medicaid: insurance through medicaid
#' - champus: military insurance
#' - hlthins: health insurance
#' - regmed: regoular source of care
#' - dri: see same doctor
#' - diabete: have diabetes
#' - hearthcond: have heart condition
#' - stroke: have stroke
#' @source JAE data archive
#' @references
#' \insertRef{KENK:TERZ:01}{micsr}
#' @importFrom Rdpack reprompt
NULL

#' Foreign exchange derivatives use by large US bank holding companies
#' @name federiv
#' @docType data
#' @keywords dataset
#' @description  a cross-section of 794 banks from 1996 to 2000 
#' @format a tibble containing:
#' - federiv: foreign exchange derivatives use, a dummy
#' - optval: option awards
#' - eqrat: leverage
#' - bonus: bonus
#' - ltass: logarithm of total assets
#' - linsown: logarithm of the percentage of the total shares outstanding that are owned by officers and directors
#' - linstown: logarithm of the percentage of the total shares outstanding that are owned by all institutional investors
#' - roe: return on equity
#' - mktbk: market to book ratio
#' - perfor: foreign to total interest income ratio
#' - dealdum: derivative dealer activity dummy
#' - div: dividends paid
#' - year: year, from 1996 to 2000
#' - no_emp: number of employees
#' - no_subs: number of subsidiaries
#' - no_off: number of offices
#' - ceo_age: CEO age
#' - gap: 12 month maturity mismatch
#' - cfa: ratio of cash flow to total assets
#' @source Lee Adkin's home page https://learneconometrics.com/
#' @references
#' \insertRef{ADKI:12}{micsr}
#' 
#' \insertRef{ADKI:CART:SIMP:07}{micsr}
#' @importFrom Rdpack reprompt
NULL

#' Demand for fees and admission from the US expense survey
#' @name feesadm
#' @docType data
#' @keywords dataset
#' @description  a cross-section of 6927 households from 2019 
#' @format a tibble containing:
#' - feesadm: expense in fees and admission
#' - linc: log of income per consumption unit
#' - smsa: a dummy indicating if the household leaves in a MSA
#' - size: number of members in the household
#' - age: the age of the head of the household
#' - educ: number of years of education of the head of the household
#' @source US bureau of labour statistics \url{https://www.bls.gov/cex/}
#' @importFrom Rdpack reprompt
NULL

#' Political economy of financial reforms
#' @name fin_reform
#' @docType data
#' @keywords dataset
#' @description  a pseudo-panel of 35 countries from 1973 to 1996 
#' @format a tibble containing:
#' - country: the country id
#' - year: the year
#' - region: the region
#' - pol: political orientation of the government
#' - fli: degree of policy liberalization index (from 0 to 18)
#' - yofc: year of office
#' - gdpg: growth rate of the gdp
#' - infl: inflation rate
#' - bop: balance of payments crises
#' - bank: banking crises
#' - imf: IMF program dummy
#' - usint: international interest rates
#' - open: trade openess
#' @source AEA website
#' @references
#' \insertRef{ABIA:MODY:05}{micsr}
#' @importFrom Rdpack reprompt
NULL

#' Food consumption in the Netherlands
#' @name food
#' @docType data
#' @keywords dataset
#' @description  a cross-section of 4611 households 
#' @format a tibble containing:
#' - year: year of the survey, either 1980 or 1988
#' - weights: weight in the survey
#' - hsize: number of persons in the household
#' - ageh: age of the head of the household in 5 classes
#' - income: total income
#' - food: food expenditure
#' - midage: a dummy for household for which the head is aged between 35 and 64 years old
#' @source Journal of Applied Econometrics Data Archive : \url{http://qed.econ.queensu.ca/jae/}
#' @references
#' \insertRef{CROM:PALM:URBA:97}{micsr}
#' @examples
#' tobit1(log(food) ~ log(income) + log(hsize) + midage,
#'        data = food, subset = year == 1980,
#'        left = -Inf, right = log(13030))
#' @importFrom Rdpack reprompt
NULL

#' Augmented Sollow's growth model with human capital
#' @name growth
#' @docType data
#' @keywords dataset
#' @description  a cross-section of 121 countries from 1985 
#' @format a tibble containing:
#' - country: country
#' - group: country's group, a factor with levels `oecd`, `other`, lqdata` (for low quality data) and `oil` (oil exporter countries)
#' - gdp60: per-capita gdp in 1960
#' - gdp85: per-capita gdp in 1985
#' - gdpgwth: mean annual growth rate of the gdp
#' - popgwth: mean annual growth rate of the population
#' - inv: mean value of the investment rate
#' - school: measurment of human capital
#' - growth: log-difference of per-capita gdp between 1985 and 1960
#' @source https://github.com/HariharanJayashankar/mrw1992
#' @references
#' \insertRef{MANK:ROME:WEIL:92}{micsr}
#' @importFrom Rdpack reprompt
NULL

#' Willingness to pay for the preservation of the Kakadu national park
#' @name kakadu
#' @docType data
#' @keywords dataset
#' @description  a cross-section of 1827 individuals from 1990 
#' @format a tibble containing:
#' - lower: lowerbound of willingness to pay, 0 if observation is left censored
#' - upper: upper bound of willingness to pay, 999 if observation is right censored
#' - recparks: the greatest value of national parks and nature reserves is in recreation activities (from 1 to 5)
#' - jobs: jobs are the most important thing in deciding how to use our natural ressources  (from 1 to 5)
#' - lowrisk: development should be allowed to proceed where environmental damage from activities such as mining is possible but very unlikely  (from 1 to 5)
#' - wildlife: it's important to have places where wildlife is preserved  (from 1 to 5)
#' - future: it's important to consider future generations  (from 1 to 5)
#' - aboriginal: in deciding how to use areas such as Kakadu national park, their importance to the local aboriginal people should be a major factor  (from 1 to 5)
#' - finben: in deciding how to use our natural ressources such as mineral deposits and forests, the most important thing is the financial benefits for Australia  (from 1 to 5)
#' - mineparks: if areas within natural parks are set aside for development projects such as mining, the value of the parks is greatly reduced  (from 1 to 5)
#' - moreparks: there should be more national parks created from state forests  (from 1 to 5)
#' - gov: the government pays little attention to the people in making decisions  (from 1 to 4)
#' - envcon: the respondent recycles things such as paper or glass and regularly buys unbleached toilet paper or environmentally friendly products ?
#' - vparks: the respondent has visited a national park or bushland recreation area in the previous 12 months ?
#' - tvenv: the respondent watchs tv programs about the environment ?  (from 1 to 9)
#' - conservation: the respondent is member of a conservation organization ?
#' - sex: male,female
#' - age: age
#' - schooling: years of schooling
#' - income: respondent's income in thousands of dollars
#' - major: the respondent received the major--impact scenario of the Kakadu conservation zone survey ?
#' @source Journal of Business, Economics and Statistics
#' @references
#' \insertRef{WERN:99}{micsr}
#' 
#' \insertRef{CARS:WILK:IMBE:94}{micsr}
#' @importFrom Rdpack reprompt
NULL

#' Resumes with 'white' or 'african-american' names
#' @name lakisha
#' @docType data
#' @keywords dataset
#' @description  a cross-section of 4870 resumes from 2001 to 2002 
#' @format a tibble containing:
#' - call: has the employer responded to the resume ?
#' - firstname: a factor of 36 first names (18 for males and 18 for females, 18 for blacks and 18 for whites)
#' - city: one of 'boston' and 'chicago'
#' - sex: one of 'female', 'male'
#' - race: one of 'black' and 'white'
#' - quality: a subjective measure of the quality of the resume, one of 'low' and 'high'
#' - exp: years of experience
#' - volunteer: volunteering experience ?
#' - military: military experience ?
#' - email: email address ?
#' - empholes: employment holes ?
#' - worksch: work in school ?
#' - honors: honors ?
#' - compsk: computer skills ?
#' - specsk: special skills ?
#' - apblack: fraction african-americans in applicant's zip code
#' - apwhite: fraction whites in applicant's zip code
#' - apdrop: fraction high school dropouts in applicant's zip code
#' - apcolp: fraction college or more in applicant's zip code
#' - alinc: log of median household income in applicant's zip code
#' - occupation: the occupation, a factor with 6 levels
#' - sector: the secor, a factor with 7 levels
#' - req: any requirement ?
#' - reqexp: experience requirement ?
#' - reqcom: communication skills required ?
#' - reqeduc: education requirement ?
#' - reqcomp: computer skills required ?
#' - reqorg: organization skills ?
#' - fed: federal contractor ?
#' - eoe: equal opportunity employer ?
#' - ownership: ownership status ?
#' - epblack: fraction african-americans in employer's zip code
#' - epwhite: fraction whites in employer's zip code
#' - epdrop: fraction high school dropouts in employer's zip code
#' - epcolp: fraction college or more in employer's zip code
#' - elinc: log of median household income in employer's zip code
#' @source American Economic Association Data Archive : \url{https://www.aeaweb.org/aer/}
#' @references
#' \insertRef{BERT:MULL:04}{micsr}
#' @importFrom Rdpack reprompt
NULL

#' US commercial loan market
#' @name loan_market
#' @docType data
#' @keywords dataset
#' @description  a time series of 72 monthly observations from 1979-01 to 1984-12 
#' @format a tibble containing:
#' - date: the month of observation
#' - loans: total commercial loans, billions of dollars
#' - prime_rate: average prime rate charged by banks
#' - aaa_rate: AAA corporate bond rate (alternative financing for firms)
#' - ipi: industrial production index (firms' expectation about future economic activity)
#' - treas_rate: 3-month treasury bill rate (alternative rate of return for banks)
#' - deposits: total bank deposit (scale variable), billions of dollars
#' @source G.S. Maddala thanks Water Mayer for providing him with the data (p.364-65)
#' @references
#' \insertRef{MADD:01}{micsr}
#' @importFrom Rdpack reprompt
NULL

#' Choice between car and transit
#' @name mode_choice
#' @docType data
#' @keywords dataset
#' @description  a cross-section of 842 individuals 
#' @format a tibble containing:
#' - mode: 1 for car, 0 for transit
#' - cost: transit fare minus automobile travel cost in US$
#' - ivtime: transit in-vehicule travel time minus in-vehicule travel time (minutes)
#' - ovtime: transit out-of vehicule time minus out-of vehicule travel time (minutes)
#' - cars: number of cars owned by the traveler's household
#' @source GAMS's website \url{http://www.gams.com/modlib/libhtml/mws.htm}
#' @references
#' \insertRef{HORO:93}{micsr}
#' @importFrom Rdpack reprompt
NULL

#'  Labor force participation and hours worked by women (National Longitudinal Survey of older women)
#' @name moffitt
#' @docType data
#' @keywords dataset
#' @description  a cross-section of 610 women from 1972 
#' @format a tibble containing:
#' - hours: hours worked per week,
#' - wage: hourly wage,
#' - nwinc: asset income per week,
#' - married: dummy for married women,
#' - age: age of the women,
#' - black: dummy for black women,
#' - fsize: size of the family
#' - clt6: number of children under 6 years old,
#' - cgt6: number of children above 6 years old,
#' - educ: number of years of education,
#' - lfsize: size of the labor force (in millions),
#' - manuf: manufacturing fraction,
#' - gov: government fraction.
#' @source this data set was kindly provided by David Drukker
#' @references
#' \insertRef{MOFF:84}{micsr}
#' 
#' \insertRef{SKEE:VELL:99}{micsr}
#' 
#' \insertRef{DRUK:02}{micsr}
#' @importFrom Rdpack reprompt
NULL

#' Mortgage Defaults
#' @name mortgage_defaults
#' @docType data
#' @keywords dataset
#' @description  yearly observations of 457 counties from 2007 to 2011 
#' @format a tibble containing:
#' - fips: county id
#' - year: year
#' - default: rate of mortgage default
#' - fico: fico
#' - impl1: impl1
#' - ur: ur
#' - loans: loans
#' @source American Economic Association Data Archive : \url{https://www.aeaweb.org/aer/}
#' @references
#' \insertRef{KUMA:18}{micsr}
#' @importFrom Rdpack reprompt
NULL

#' National supported work demonstration
#' @name nsw
#' @docType data
#' @keywords dataset
#' @description  a cross-section of 19204 individuals from 1978 
#' @format a tibble containing:
#' - smpl: the sample, one of `nws` (the experimental NWS sample), `psid` (PSID control group) and `cps` (CPS control group)
#' - sub: an integer to select a subsample used by Dehejia and Wahba: `sub == 1` select the whole sample, `smpl == "nws"
#' - group: one of `"control"` and `"treated"`
#' - ethnicity: a factor with levels `"other"`, `"black"` and `"hispanic"`
#' - age: ag in years
#' - education: years of education
#' - married: a dummy for married individuals
#' - nodegree: a dummy for individuals who left high school without a degre
#' - re74: income in 1974
#' - re75: income in 1975
#' - re78: income in 1978
#' @source `https://users.nber.org/~rdehejia/data/.nswdata2.html`
#' @references
#' \insertRef{LALO:86}{micsr}
#' 
#' \insertRef{DEHE:WAHB:99}{micsr}
#' 
#' \insertRef{DEHE:WAHB:02}{micsr}
#' @importFrom Rdpack reprompt
NULL

#' Oil investment
#' @name oil
#' @docType data
#' @keywords dataset
#' @description  a cross-section of 53 oil fields from 1969 to 1992 
#' @format a tibble containing:
#' - dur: duration of the appraisal lag in mounths (time span between discovery of an oil field and beginning of development, i.e. approval of annex B).
#' - size: size of recoverable reserves in millions of barrels
#' - waterd: depth of the sea in metres
#' - gasres: size of recoverable gas reserves in billions of cubic feet
#' - operator: equity market value (in 1991 million pounds) of the company operating the oil field
#' - p: real after--tax oil price measured at time of annex B approval
#' - vardp: volatility of the real oil price process measured as the squared recursive standard errors of the regression of pt-pt-1 on a constant
#' - p97: adaptative expectations (with parameter theta=0.97) for the real after--tax oil prices formed at the time of annex B approval
#' - varp97: volatility of the adaptative expectations (with parameter theta=0.97) for real after tax oil prices measured as the squared recursive standard errors of the regression of pt on pte(theta)
#' - p98: adaptative expectations (with parameter theta=0.98) for the real after--tax oil prices formed at the time of annex B approval
#' - varp98: volatility of the adaptative expectations (with parameter theta=0.98) for real after tax oil prices measured as the squared recursive standard errors of the regression of pt on pte(theta)
#' @source Journal of Applied Econometrics Data Archive : \url{http://qed.econ.queensu.ca/jae/}
#' @references
#' \insertRef{FAVE:PESA:SHAR:94}{micsr}
#' @importFrom Rdpack reprompt
NULL

#' Voucher and private school
#' @name paces
#' @docType data
#' @keywords dataset
#' @description  a cross-section of 1577 households from 1995 
#' @format a tibble containing:
#' - id: NA
#' - privsch: NA
#' - educyrs: NA
#' - voucher: NA
#' - pilot: NA
#' - housvisit: NA
#' - city: NA
#' - phone: NA
#' - age: NA
#' - sex: NA
#' - strata: NA
#' - smpl: NA
#' - month: NA
#' - married: NA
#' - finish8: V
#' - repetitions: NA
#' - in_school: NA
#' - year: NA
#' @source Joshua Angrist's web site
#' @references
#' \insertRef{ANGR:BETT:BLOO:KING:KREM:02}{micsr}
#' @importFrom Rdpack reprompt
NULL

#' Precautionary motives and portfolio decisions
#' @name portfolio
#' @docType data
#' @keywords dataset
#' @description  yearly observations of 3348 households from 1993 to 1998 
#' @format a tibble containing:
#' - id: household's identifier
#' - year: year
#' - share: share of riskless assets
#' - uncert: uncertainty, one of low, moderate, high,
#' - expinc: income expectation, one of increase, constant and decrease
#' - finass_10: log of financial assets < 10k Dfl.
#' - finass_10_100: log of financial assets 10 - 100k Dfl.
#' - finass_more: log of financial assets > 100k Dfl.
#' - networth: log of net worth
#' - noncapinc: log of non-capital income
#' - mtrate: Household marginal tax rate
#' - high_inc_oversmpl: High-income oversample
#' - age: age
#' - educ: one of low, second, vocat and univ
#' - diploma: dummy
#' - female: dummy
#' - adults: number of adults
#' - child_0_12: number of children under 13 years old
#' - child_13_more: number of children of 13 years or more
#' - occup: occupation, one of paidjob, unempl, retired, disab, self and other
#' - riskav: degree of risk aversion, one of high, low, medium and dontknow
#' - feeling: one of veryhappy, happy, neither (happy or unhappy), unhappy
#' - flex: job flexibility, one of none, addjob, hours and both
#' - smoke: smoking habits, one of none, low and high
#' - alcohol: daily alcohol, dummy for a consumption of at least 4 drinks
#' - body_mass: body mass index
#' - habits: one of none, strong and weak
#' @source Journal of Applied Econometrics Data Archive : \url{http://qed.econ.queensu.ca/jae/}
#' @references
#' \insertRef{HOCH:03}{micsr}
#' @importFrom Rdpack reprompt
NULL

#' Market shares of train and air for Paris and main French cities in 1995
#' @name price_time
#' @docType data
#' @keywords dataset
#' @description  a cross-section of 13 cities from 1995 
#' @format a tibble containing:
#' - town: the city
#' - region: one of `SW`, `SE` and `O`
#' - trafic_rail: number of rail passengers
#' - trafic_air: number of air passengers
#' - price_rail: price of a train ticket in euros
#' - price_air: price of a plane ticket in euros
#' - time_rail: length of a train trip in minutes
#' - time_air: length a plane trip in minutes
#' @source Patrick Bonnel's book
#' @references
#' \insertRef{BONN:04}{micsr}
#' @importFrom Rdpack reprompt
NULL

#' Private insurance and doctor visits
#' @name private_ins
#' @docType data
#' @keywords dataset
#' @description  a cross-section of 7555 individuals from 2010 
#' @format a tibble containing:
#' - size: family size
#' - smsa: lives in a msa
#' - age: age
#' - sex: sex
#' - educ: number of years of education
#' - wage: wage
#' - nbemp: number of employees in the firm
#' - pluriloc: a dummy for firms with multiple locations
#' - privateins: a dummy for private insurance
#' - doctor: a dummy if doctor visits
#' - married: a dummy equal 1 if ever married
#' - income: income
#' - region: a factor for the region, one of `northeast`, `midwest`, `south` and `west`
#' - race: a factor for the race, one of `white`, `black`, `minority` and `asian`
#' - mental: a factor for mental health, one of `poor`, `fair`, `good`, `verygood` and `excellent`
#' - physical: a factor for physical health, one of `poor`, `fair`, `good`, `verygood` and `excellent`
#' @source JAE data archive
#' @references
#' \insertRef{HAN:LEE:19}{micsr}
#' @importFrom Rdpack reprompt
NULL

#' Canadian Manufacturing
#' @name prod_canada
#' @docType data
#' @keywords dataset
#' @description  a time series of 25 yearly observations from 1946 to 1970 
#' @format a tibble containing:
#' - year: year
#' - po: output price
#' - qo: output volume
#' - pm: price index for materials
#' - qm: quantity of materials
#' - pl: price index of productive labor
#' - ql: quantity of productive labor
#' - pn: price index of non-productive labor
#' - qn: quantity of non-productive labor
#' - k: capital stock
#' @source original paper, p.375
#' @references
#' \insertRef{WOOD:77}{micsr}
#' @importFrom Rdpack reprompt
NULL

#' Education of oldest sibling as an instrument of own education
#' @name sibling_educ
#' @docType data
#' @keywords dataset
#' @description  a cross-section of 303 individuals from 1992 
#' @format a tibble containing:
#' - wage: hourly wage
#' - educ: Education (years of schooling)
#' - ability: ability based on standardized AFQT test score
#' - educsib: education of oldest sibling (years of schooling)
#' - experience: experience (total weeks of labor market experience)
#' - tenure: tenure (weeks on current job)
#' - mothed: Mother's education (years of schooling)
#' - fathed: Father's education (years of schooling)
#' - residence: one of urban or rural
#' - broken: broken home dummy (dummy for living with both parents at age 14)
#' @source Journal of Applied Econometrics Data Archive : \url{http://qed.econ.queensu.ca/jae/}
#' @references
#' \insertRef{KOOP:POIR:TOBI:05}{micsr}
#' @importFrom Rdpack reprompt
NULL

#' Number of direct job changes from the German Socio-economic panel survey
#' @name soep
#' @docType data
#' @keywords dataset
#' @description  a cross-section of 2651 individuals from 1984 
#' @format a tibble containing:
#' - djc: number of direct job changes, i.e. changes without an intervening spell of unemployment,
#' - lfp: a dummy for labor force participation,
#' - exper: full-time employment in years
#' - single: a factor indicating if the individual is single
#' - whitecol: a foctor indicating if the individual is white collar
#' - educ: lenght of education in years
#' - spd: a factor indicating if the individual is a strong or a very strong supporter of the social-democrat party of Germany
#' @source this data is part of `SemiParSampleSel` package
#' @references
#' \insertRef{WYSZ:MARR:18}{micsr}
#' @importFrom Rdpack reprompt
NULL

#' Solow's growth model with spatial correlation
#' @name sp_solow
#' @docType data
#' @keywords dataset
#' @description  a cross-section of 91 countries from 1995 
#' @format a tibble containing:
#' - name: the name of the country
#' - code: the id of the country
#' - capital: name of the capital
#' - gdp60: per capita gdp in 1960
#' - gdp95: per capita gdp in 1995
#' - saving: saving rate
#' - labgwth: growth rate of the labor force
#' - point: the coordinates of the capital
#' - border: the coordinates of the country
#' @source JAE data archive
#' @references
#' \insertRef{ERTU:KOCH:07}{micsr}
#' @importFrom Rdpack reprompt
NULL

#' Railroad tracks and segregation
#' @name tracks_side
#' @docType data
#' @keywords dataset
#' @description  a cross-section of 121 cities from 1990 
#' @format a tibble containing:
#' - city: city
#' - state: state
#' - segregation: dissimilarity index in 1990 of the city
#' - trackdens: track length per square kilometer
#' - raildiv: railroad division index
#' - giniw: gini index for whites
#' - ginib: gini index for blacks
#' - povw: poverty rate for whites
#' - povb: poverty rate for blacks
#' - D9wD9b: ratio of the ninth decile for whites and blacks
#' - D1wD1b: ratio of the first decile for whites and blacks
#' - D9wD1b: ratio of the ninth decile for whites and the first decile for blacks
#' - D9bD1w: ratio of the ninth decile for blacks and the first decile for whites
#' @source American Economic Association Data Archive : \url{https://www.aeaweb.org/aer/}
#' @references
#' \insertRef{ANAN:11}{micsr}
#' @importFrom Rdpack reprompt
NULL

#' Lobying from Capitalists and Unions and Trade Protection
#' @name trade_protection
#' @docType data
#' @keywords dataset
#' @description  a cross-section of 194 United States 
#' @format a tibble containing:
#' - ntb: NTB coverage ratio, proportion
#' - exports: exportations
#' - imports: importations
#' - elast: demand elasticity
#' - cap: lobying
#' - labvar: labor market covariate
#' - sic3: 3-digit SIC industry classification
#' - k_serv: bla
#' - inv: Inventories, factor share
#' - engsci: Engineers and scientists, factor share
#' - whitecol: White collar, factor share
#' - skill: Skilled, factor share
#' - semskill: Semi-skilled, factor share
#' - cropland: Cropland, factor shaer
#' - pasture: Pasture, factor share
#' - forest: Forest, factor share
#' - coal: Coal, factor share
#' - petro: Petroleum, factor share
#' - minerals: Minerals, factor share
#' - scrconc: Seller concentration
#' - bcrconc: Buyer concentration
#' - scrcomp: Seller number of firms
#' - bcrcomp: Buyer number of firms
#' - meps: Scale
#' - kstock: Capital stock
#' - puni: bla
#' - geog2: Geographic concentration
#' - tenure: Average worker tenure, years
#' - klratio: Capital-labor ratio
#' - bunion: bla
#' @source American Economic Association Data Archive : \url{https://www.aeaweb.org/aer/}
#' @references
#' \insertRef{MATS:SHER:06}{micsr}
#' @importFrom Rdpack reprompt
NULL

#' Political Economy and Traffic Citations
#' @name traffic_citations
#' @docType data
#' @keywords dataset
#' @description  a cross-section of 68357 fines from 2001 
#' @format a tibble containing:
#' - fine: \code{yes} for a fine, \code{no} for a warning,
#' - amount: the amount of the fine if issued
#' - black: dummy for a black offender,
#' - hispanic: dummy for an hispanic offender,
#' - age: the age of the offender
#' - female: dummy for a female offender,
#' - cdl: dummy for a commercial driving license,
#' - res: the residence of the offender, one of \code{loc}, \code{oto} and \code{ost}, whether the offender resides in the town, out of the town in Massaschusetts or in another state
#' - courtdist: distance from the distance courthouse
#' - id: the id of the police officer
#' - stpol: dummy for state officers
#' - mph: the difference between actual and limit speed,
#' - prval: the average property value of the municipality
#' - hospemp: the percentage of employement in the hospitality sector
#' - oloss: a dummy that indicates that the municipality failed to pass an overide referendum
#' @source American Economic Association Data Archive : \url{https://www.aeaweb.org/aer/}
#' @references
#' \insertRef{MAKO:STRA:09}{micsr}
#' @examples
#' sel <- fine ~ log(mph) + black + hispanic + female * log(age) +
#'        stpol * (res + log(prval) + oloss) + cdl
#' out <- update(sel, log(amount) ~ . - cdl)
#' probit <- glm(sel, data = traffic_citations,
#'              family = binomial(link='probit'))
#' lp <- probit$linear.predictor
#' mls <- dnorm(lp) / pnorm(lp)
#' summary(probit)
#' # Computing the heckit by hand
#' lm <- lm(out, traffic_citations)
#' heck <- update(lm, . ~ . + mls)
#' # or using the heckit function
#' library("sampleSelection")
#' heck <- heckit(sel, out, data = traffic_citations,
#'             method = "2step")
#' summary(heck)
#' @importFrom Rdpack reprompt
NULL

#' Determinants of household trip taking
#' @name trips
#' @docType data
#' @keywords dataset
#' @description  a cross-section of 577 households from 1978 
#' @format a tibble containing:
#' - trips: number of trips taken by a member of a household the day prior the survey interview
#' - car: 1 if household owns at least one motorized vehicule
#' - workschl: share of trips for work or school vs personal business or pleasure
#' - size: number of individuals in the household
#' - dist: distance to central business district in kilometers
#' - smsa: a factor with levels "small" (less than 2.5 million population) and "large" (more than 2.5 million population)
#' - fulltime: number of fulltime workers in household
#' - adults: number of adults in household
#' - distnod: distace from home to nearest transit node, in blocks
#' - realinc: household income divided by median income of census tract in which household resides
#' - weekend: 1 if the survey period is either saturday or sunday
#' @source kindly provided by Joseph Terza
#' @references
#' \insertRef{TERZ:98}{micsr}
#' 
#' \insertRef{TERZ:WILS:90}{micsr}
#' @importFrom Rdpack reprompt
NULL

#' Temporary help jobs and permanent employment
#' @name twa
#' @docType data
#' @keywords dataset
#' @description  a cross-section of 2030 individuals 
#' @format a tibble containing:
#' - id: identification code
#' - age: age
#' - sex: a factor with levels `"female"` and `"male"`
#' - marital: marital status, `"married"` or `"single"`
#' - children: number of children
#' - feduc: father's education
#' - fbluecol: father blue-color
#' - femp: father employed at time 1
#' - educ: years of education
#' - pvoto: mark in last degree as fraction of max mark
#' - training: received professional training before treatment
#' - dist: distance from nearest agency
#' - nyu: fraction of school-to-work without employment
#' - hour: weekly hours of work
#' - wage: monthly wage
#' - hwage: hourly wage at time 1
#' - contact: contacted a temporary work agency
#' - region: one of  `"Tuscany"` and `"Sicily"`
#' - city: the city
#' - group: one of `"control"` and `"treated"
#' - sector: NA
#' - occup: occupation, one of `"nojob"`, `"selfemp"`, `"bluecol"` and `"whitecol"`
#' - empstat: employment status, one of `"empl"`, `"unemp"` and `"olf"` (out of labor force)
#' - contract: job contract, one of `"nojob"`, `"atyp"` (atypical) and `"perm"` (permanent)
#' - loc: localisation, one of `"nord"`, `"centro"`, `"sud"` and `"estero"`
#' - outcome: one of `"none"`, `"other"`, `"fterm"` and `"perm"`
#' @source Journal of Applied Econometrics Data Archive : \url{http://qed.econ.queensu.ca/jae/}
#' @references
#' \insertRef{ICHI:MEAL:NANN:08}{micsr}
#' @importFrom Rdpack reprompt
NULL

#' Education and earnings of british twins
#' @name twins
#' @docType data
#' @keywords dataset
#' @description  a pseudo-panel of 214 families from 1999 
#' @format a tibble containing:
#' - family: family index
#' - earning: hourly wage rate
#' - educ: reported years of schooling
#' - educt: estimated years of schooling
#' - age: age
#' - londse: lives in
#' - part: part-time job ?
#' - married: married ?
#' - tenure: tenure in years
#' @source American Economic Review data archive : \url{http://aeaweb.org/aer/}
#' @references
#' \insertRef{BONJ:CHER:HASK:HAWK:SPEC:03}{micsr}
#' @importFrom Rdpack reprompt
NULL

#' Density and distance to the CBD in Alabama
#' @name urban_gradient
#' @docType data
#' @keywords dataset
#' @description  a cross-section of 2315 blocks from 1999 
#' @format a tibble containing:
#' - msa: msa name
#' - county: county name
#' - tract: tract id
#' - blkg: blkg id
#' - area: area in square km
#' - population: number of inhabitants
#' - distance: distance to the CBD in km
#' @source American Economic Association Data Archive : \url{https://www.aeaweb.org/aer/}
#' @references
#' \insertRef{DURA:PUGA:20}{micsr}
#' @importFrom Rdpack reprompt
NULL

#' US Counties
#' @name us_counties
#' @docType data
#' @keywords dataset
#' @description  a cross-section of 3141 counties 
#' @format a tibble containing:
#' - fips: fips code
#' - gid: gadm code
#' - state: state name
#' - county: county name
#' - geometry: a sfc containing the coordinates of counties border
#' @source gadm of the geodata package
#' @importFrom Rdpack reprompt
NULL

#' Draft eligibility and income
#' @name vn_income
#' @docType data
#' @keywords dataset
#' @description  a cross-section of 40880 cells of individuals 
#' @format a tibble containing:
#' - birth: year of birth
#' - year: year of observation for the income
#' - interval: 5 days interval of date of birth (from 1 to 73)
#' - type: income measurment, one of `fica` and `w2`
#' - race: a factor with levels `white` and `nonwhite`
#' - freq: the number of individuals in the cell
#' - null_inc: share of null incomes in the cell
#' - mean_inc: mean income in the cell
#' - sd_inc: standard-deviation of income in the cell
#' @source Joshua Angrist's web page [](https://economics.mit.edu/faculty/angrist)
#' @references
#' \insertRef{ANGR:90}{micsr}
#' @importFrom Rdpack reprompt
NULL

#' Draft eligibility and enrolment in the Viet-Nam war
#' @name vn_veteran
#' @docType data
#' @keywords dataset
#' @description  a cross-section of 4545 individuals 
#' @format a tibble containing:
#' - birth: year of birth
#' - wgts: weights
#' - race: a factor with levels `white` and `nonwhite`
#' - veteran: a factor with levels `yes` for veterans and `no` otherwise
#' - eligible: a factor with levels `yes` for draft eligible men and `no` otherwise
#' @source Joshua Angrist's web page [](https://economics.mit.edu/faculty/angrist)
#' @references
#' \insertRef{ANGR:90}{micsr}
#' @importFrom Rdpack reprompt
NULL

#' Price of wine
#' @name wine
#' @docType data
#' @keywords dataset
#' @description  a cross-section of 9600 wines 
#' @format a tibble containing:
#' - price: wine price (cpi adjusted)
#' - cases: cases produced
#' - score: wsm tasting score
#' - age: years of aging
#' - region: region of production
#' - variety: grape variety used
#' - vintage: vintage, from 1991 to 1999
#' - reserve: \code{'yes'} if reserve wine
#' - vineyard: \code{'yes'} if vineyard information is provided
#' - estate: \code{'yes'} if estate produced
#' - class: cluster 1,2,3,4, per LPRC
#' @source Journal of Applied Econometrics Data Archive : \url{http://qed.econ.queensu.ca/jae/}
#' @references
#' \insertRef{COST:MITT:MCCL:09}{micsr}
#' @importFrom Rdpack reprompt
NULL

#' Economies and diseconomies of scale for US counties
#' @name agglo_growth
#' @docType data
#' @keywords dataset
#' @description  a cross-section of 3106 counties from 1990 
#' @format a tibble containing:
#' - fips: county fips code
#' - emp_gr: employment growth rate
#' - pop_gr: population growth rate
#' - emp: employment in 1980
#' - pop: population in 1980
#' - college: fraction of adult population with bachelor's degree or more in 1980
#' - manuf: fraction of employment in manufacturing in 1980
#' - unemp: unemployment rate in 1980
#' - income: per capita income in 1979
#' - educ_sh: share of local government spending on education in 1982
#' - hw_sh: share of local government spending on highways in 1982
#' - pol_sh: share of local government spending on police in 1982
#' - notwhite: fraction of population that is not white in 1980
#' - type: a factor with levels `"rural"` and `"urban"`
#' @source JAE data archive
#' @references
#' \insertRef{WHEE:03}{micsr}
#' @importFrom Rdpack reprompt
NULL
