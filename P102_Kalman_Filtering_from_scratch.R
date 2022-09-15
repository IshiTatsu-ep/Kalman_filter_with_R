#P.102　カルマンフィルタリング

#ナイル川の流域データを使用
y <- Nile
t_max <- length(y)

#1時点分のカルマンフィルタリングを行う関数
Kalman_filtering <- function(m_t_minus_1,C_t_minus_1,t){
  #一期先予測分布
  a_t <- G_t %*% m_t_minus_1
  R_t <- G_t %*% C_t_minus_1 %*% t(G_t) + W_t
  
  #一期先予測尤度
  f_t <- F_t %*% a_t
  Q_t <- F_t %*% R_t %*% t(F_t) +V_t
  
  #カルマンゲイン
  K_t <- R_t %*% t(F_t) %*% solve(Q_t)
  
  #状態の更新
  m_t <- a_t + K_t %*% (y[t]-f_t)
  C_t <- diag(nrow(R_t) - K_t %*% F_t) %*% R_t
  
  #フィルタリング分布と一期先予測分布の平均と分散を返す
  return(list(m = m_t, C= C_t, 
              a = a_t, R = R_t))
}

#線形ガウス型状態空間のパラメータを設定（すべて1×1の行列）
G_t <- matrix(1, ncol=1, nrow=1); W_t <- matrix(exp(7.29), ncol=1, nrow=1)
F_t <- matrix(1, ncol=1, nrow=1); V_t <- matrix(exp(9.62), ncol=1,nrow=1)
m0 <- matrix(0,ncol=1,nrow=1); C0 <- matrix( 1e+7, ncol=1,nrow=1)

#フィルタリング分布（+一期先予測分布）の平均分散を求める
#状態ベクトルの初期化
m <- rep(NA_real_, t_max) ; C <- rep(NA_real_, t_max)
a <- rep(NA_real_, t_max) ; R <- rep(NA_real_, t_max)

#時点：t-1
KF <- Kalman_filtering(m0,C0,t = 1)
m[1] <- KF$m; C[1] <- KF$C
a[1] <- KF$a; R[1] <- KF$R

#時点2～t_max
for(t in 2:t_max){
  KF <- Kalman_filtering(m[t-1],C[t-1],t=t)
  m[t] <- KF$m; C[t] <- KF$C
  a[t] <- KF$a; R[t] <- KF$R
}

y_plot <- as.vector(y)
#out <- cbind(y_plot,m)
#ts.plot(y_plot)
ts.plot(cbind(y_plot,m),
        lty = c("solid","solid"),
        col = c("lightgray", "blue"),
        ylab = "")

#凡例
legend(legend = c("観測値","カルマン平滑化"),
       lty = c("solid","solid"),
       col = c("lightgray","blue"),
       x = "topright", text.width = 17, cex=0.6)




