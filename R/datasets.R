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
#' - dindx: difference of the inflation rate
#' - indx: inflation rate divided by 18
#' - indxl: lag value of indx
#' - rhs1: indxl * (1 - indxl)
#' - max_indxl: maximumum value of indxl by year and region
#' - catchup: difference between max_indxl and indxl
#' - dum_bop: balance of paiement crisis in the first two previous years
#' - dum_bank: bank crises in the first two previous years
#' - dum_1yofc: dummy for first year of office
#' - recession: dummy for recessions
#' - hinfl: dummy for inflation rate greater than 50 percent
#' @source AEA website
#' @references
#' \insertRef{ABIA:MODY:05}{micsr}
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
#' @source GAMS's website \url{https://www.gams.com/latest/gamslib_ml/libhtml/gamslib_mws.html}
#' @references
#' \insertRef{HORO:93}{micsr}
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
#' - sector: the sector
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

