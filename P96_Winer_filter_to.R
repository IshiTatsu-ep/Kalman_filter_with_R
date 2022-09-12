#P.96 AR(1)モデルにおけるウィーナー平滑化とカルマン平滑化

#前処理
set.seed(23)
#install.packages('dlm',dependencies = TRUE)
library(dlm)

#AR(1)を含む状態空間の設定
W <- 1
V <- 2
phi <- 0.98
#dlmModPolyは任意のorderを持つ多項式モデルを設定できる。
#dw = システムノイズの分散行列
#dv = 観測ノイズの分散
#C0 = 状態ベクトルのサンプリング前の分散行列
mod <- dlmModPoly(order = 1, dW = W, dV = V, C0 = 100)
mod$GG[1,1] <- phi

#カルマン予測を活用して、観測値（ダミーデータ）を作成
t_max <- 100
sim_data <- dlmForecast(mod = mod, nAhead = t_max, sampleNew = 1)
y<-sim_data$newObs[[1]]

#カルマン平滑化
dlmSmoothed_obj <- dlmSmooth(y = y, mod = mod)
s <- dropFirst(dlmSmoothed_obj$s)

#ウィナー平滑化
#係数の設定
r <- V/W
b <- 1/(r*phi) + 1/phi + phi
beta <- (b - sqrt(b^2 - 4))/2

#観測値が有限のため、前後に必要最低限のダミー0を補填
y_expand <- c(rep(0,t_max-1),y,rep(0,t_max-1))

#ウィナー平滑化の実行
d <- (1/phi - beta)*(phi - beta)/(1-beta^2) * 
  filter(method = "convolution", filter = beta^abs(-(t_max-1):(t_max-1)),
         x = y_expand)

#dからダミー成分のNAを除去
d <- d[!is.na(d)]

ts.plot(cbind(y,d,s),
        lty = c("solid","dashed","solid"),
        col = c("lightgray", "red", "blue"),
        ylab = "")

#凡例
legend(legend = c("観測値", "ウィナー平滑化", "カルマン平滑化"),
       lty = c("solid","dashed","solid"),
       col = c("lightgray", "red", "blue"),
       x = "topright", text.width = 17, cex=0.6)