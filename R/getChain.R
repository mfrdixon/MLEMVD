#' Load the option chain from a csv file into a dataframe
#'
#' @param fileName Filename of the csv file storing
#' @export output A dataframe representing the znga option chain


getChain<-function(fileName){

  ZNGA_f<-load(fileName)

  keys<-unique(ZNGA_f$TIMESTAMP)
  n<- dim(ZNGA_f)[1]

  methods::setClass("Option", methods::representation(symbol = "character", timestamp="character", underlying = "numeric", strike = "numeric", mid_price = "numeric", maturity = "numeric", type="character"))

  new.hashtable <- function() {
    e <- new.env()
    list(set = function(key, value) assign(as.character(key), value, e),
         get = function(key) get(as.character(key), e),
         rm = function(key) rm(as.character(key), e))
  }
  ht <- new.hashtable()

  for (i in 1:n){
    key <- as.character(ZNGA_f$TIMESTAMP[i])
    symbol <- as.character(ZNGA_f$SYMBOL[i])
    strike <- ZNGA_f$STRIKE[i]
    underlying <- ZNGA_f$UNDERLYING[i]
    mid_price <- (ZNGA_f$BID[i] + ZNGA_f$ASK[i])/2.0
    type  <- as.character(ZNGA_f$TYPE[i])
    t<-as.POSIXct(as.numeric(as.character(ZNGA_f$TIMESTAMP[i])), origin='1970-01-01', tz='GMT')
    T<-as.POSIXct(as.numeric(as.character(ZNGA_f$MATURITY[i])), origin='1970-01-01', tz='GMT')

    maturity <- as.integer(T-t)/252

    if (is.na(T-t)==FALSE) {
      option<- new("Option", symbol = symbol,
          underlying = underlying,
          timestamp = key,
          strike = strike,
          mid_price = mid_price,
          maturity = maturity,
          type = type)
      tryCatch(ht$get(key), error = function(e) {
           ht$set(key, list(option));
        }, finally = {chain<-ht$get(key);
                      chain<-append(chain,option);
                      ht$set(key,chain);
        })
    }
  }

  SM<-sort(unique(ZNGA_f$MATURITY))[1]

  atm_list<-list()
  for (key in keys){
   min_dist<-999
   chain <- ht$get(key)
   t<-as.POSIXct(as.numeric(key), origin='1970-01-01', tz='GMT')
   T<-as.POSIXct(as.numeric(SM), origin='1970-01-01', tz='GMT')
   maturity <- as.integer(T-t)/252
   atm_el <- NULL
   for (el in chain){
    dist<-abs(el@strike -el@underlying)

    if (dist<min_dist && el@maturity <=maturity && el@type=='C'){
      min_dist = dist;
      atm_el = el;
    }
   }
   atm_list<-append(atm_list, atm_el)
  }

  df<-data.frame(matrix(ncol=7, nrow=length(keys)))
  names(df)<-c('timestamp', 'symbol', 'underlying', 'type', 'maturity', 'strike', 'mid_price')
  i<-1
  for (el in atm_list)
  {
    df[i,'timestamp']  <- el@timestamp
    df[i,'symbol']     <- el@symbol
    df[i,'underlying'] <- el@underlying
    df[i,'type']       <- el@type
    df[i,'maturity']   <- el@maturity
    df[i,'strike']     <- el@strike
    df[i,'mid_price']  <- el@mid_price
    i <- i+1
  }
return(df)
}
