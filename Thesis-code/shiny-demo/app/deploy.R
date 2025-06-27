library(rsconnect)
rsconnect::setAccountInfo(name='msgplin',
                          token='3FA378B4F910C81AF757164F7E1B10BD',
                          secret='9WhrnviGqkv47XB9oL85tqrVG8CvpemrTjBVLSvw')
rsconnect::deployApp(
  appDir = "~/Rcode/Paper.code/app",
  appName = "Model-Discrimination-Design"
)