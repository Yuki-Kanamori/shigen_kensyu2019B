## 2019年資源管理研修（初級）
## 2019年11月19日（火）
## MSYと新しいABC算定ルール

## 簡単なシミュレーションプログラムを実行することにより新しいMSYの考え方やABCルールのポイントを理解する

## MSYの計算

# packageの読み込み
library(tidyverse)    # データ成形やグラフ作成など多岐に渡るデータを扱うパッケージを統合したもの
library(gridExtra)    # ggplotで作った異なる図を一緒に表示するパッケージ
library(VGAM)    # 今回の研修では，特に必要ないが，SRfunctionsで必要になる．

# ここは，加入変動を考えたMSYがどういう性質を持っていて，今までの（古典的な）MSYとどのように違うかを簡単なシミュレーションで見てみよう，という意図

a <- 1.2    # 再生産曲線の増加率
b <- 1000   # 再生産曲線の密度効果
Y <- 500    # MSY計算の年数
M <- 0.4    # 自然死亡係数
S <- exp(-M)   # 生残率

SR <- "HS"   # 再生産曲線がHockey-Stickか（HS）かBeverton-Holt （BH）か選択

# 再生産曲線
HS <- function(x,a=1.2,b=1000) pmin(a*x,a*b)
BH <- function(x,a=1.2,b=1000) a*x/(1+x/b)

# HSとBHを引数で選べるようにした再生産曲線
SRfunc <- function(x,a=1.2,b=1000,SR="HS"){
  if (SR=="HS") output <- HS(x,a,b)
  if (SR=="BH") output <- a*x/(1+x/b)
  
  output
}

if (SR=="HS"){ # 再生産曲線がHSの場合の決定論的FmsyやN0
  Fmsy <- 1-1/(a+S)
  Nmsy <- b/(1-Fmsy)
  R0 <- a*b
  N0 <- R0/(1-S)
}
if (SR=="BH"){ # 再生産曲線がBHの場合の決定論的FmsyやN0
  Fmsy.Est <- function(x) (b/(1-x) + b*x/(1-x)^2)*(a*(1-x)/(1-S*(1-x))-1)-b*x/(1-x)*(a/(1-S*(1-x))+a*(1-x)*S/(1-S*(1-x))^2)    #  BHは微分すると複雑な式になる
  Fmsy <- uniroot(Fmsy.Est,c(0,1))$root   # 上の関数（dC/dF）が0になるものがFmsy
  Nmsy.Est <- function(x) b/(1-x)*(a*(1-x)/(1-S*(1-x))-1)    # 平衡状態の資源量の式
  Nmsy <- Nmsy.Est(Fmsy)
  N0 <- Nmsy.Est(0)    # 平衡状態の資源量の式にF=0を与えれば，漁獲開始前の資源量になる
  R0 <- N0*(1-S)   # 漁獲なしの平衡状態N=SN+Rから，R=(1-S)Nとなる
}
MSY <- Fmsy*Nmsy

# 基本的に，SR="HS"で計算を行うが，SR="BH"とすることにより，真の再生産関係がBHのときに推定モデルとしてHSを使ったら，どれだけインパクトがあるかを見ることができる．

# MSYを計算するためのシンプルなシミュレーション
Sim.msy <- function(Fr=0, Y=500, sigma=0.2, a=1.2, b=1000, seed=1){
  # Fr：漁獲率
  # Y：シミュレーションの年数
  # sigma：加入変動の大きさ
  # a：再生産関係の増加率パラメータ
  # b：再生産関係の密度効果パラメータ
  # seed：乱数の設定
  # シミュレーションをする場合，通常これを何回か繰り返すが，簡単のため1回だけ実行するようにしている
  
  set.seed(seed)
  N <- C <- numeric(Y)
  N[1] <- Nmsy
  z <- rnorm(Y)    # 加入変動のもとになるランダム誤差
  
  # 個体群動態モデルはN_{t+1}=S*N_t*(1-Fr) + R_{t+1}という形．R_{t+1} = HS(N_t*(1-Fr))*errorとなる．漁獲された後の資源量が再生産に寄与するという仮定を用いている．この場合，N*(1-F)が親魚量ということになる．より一般的なモデルでは，成長の効果も入れるが，簡単のため成長を考えていない．
  for (i in 1:Y){
    C[i] <- N[i]*Fr     # 漁獲量の計算
    N[i+1] <- pmax(S*(N[i]-C[i])+SRfunc(N[i]-C[i],a=a,b=b,SR=SR)*exp(sigma*z[i]-sigma^2/2),0)   # 対数正規誤差を仮定（バイアスを補正）
  }
  
  list(C=C,N=N)    # 漁獲量と資源量をアウトプットする
}

# 上のシミュレーションプログラムを使って，加入の変動の大きさをいろいろ変えた場合の資源量や漁獲量を計算する
Res.sMSY <- list()    # シミュレーションの結果をとっておくために入れ物を作っておく
NY.curve <- NULL    # シミュレーション結果を要約したものをとっておく入れ物
k <- 1
Fr.range <- seq(0,0.99,by=0.01)    # 漁獲率を0から0.99まで変えて，上のシミュレーションを実行する．最大漁獲量が得られるFrはFmsyである．本来は，平衡状態に達したもので考えるが，簡単のためその辺の細かいことは無視している．
for (CV in c(0,0.2,0.4,0.6,0.8,1.0)){    # 加入の変動係数を0（決定論的なモデル），0.2，0.4，0.6，0.8，1と変えてやってMSYの変化を見る
  Res.sMSY[[k]] <- map(Fr.range, Sim.msy, sigma=CV)   # Fを変えたシミュレーション結果をRes.sMSYに格納
  NY.curve <- rbind(NY.curve, cbind(CV,Fr.range,map_dbl(1:100,function(i) mean(Res.sMSY[[k]][[i]]$N)),map_dbl(1:100,function(i) mean(Res.sMSY[[k]][[i]]$C))))    # 上の結果から，Fを変えたときの平均資源量，平均漁獲量を記録しておく
  k <- k + 1
}

NY.curve <- as.data.frame(NY.curve)    # data.frameにしてやってデータの取り扱いを容易にする
colnames(NY.curve) <- c("CV","F","Abundance","Yield")    # データの列にラベルをつけてやる
NY.curve[,1] <- factor(NY.curve[,1])    # 1列目のCVは数字だが，factorとしてカテゴリ変数だとしてやる

# MSYの計算
MSY.stat <- NY.curve %>%
  group_by(CV) %>%
  summarise(Fmsy=F[which.max(Yield)],MSY=max(Yield),Nmsy=Abundance[which.max(Yield)])    # CVごとに，Fmsy，MSY，Nmsyを計算してまとめた表にしてやる

# PGY (Pretty Good Yield）を計算
PGY.stat <- (NY.curve %>% group_by(CV) %>% filter(Yield >= 0.8*max(Yield))) %>% summarize(Flow=F[Yield==first(Yield)],PGYlow=first(Yield),Fhigh=F[Yield==last(Yield)],PGYhigh=last(Yield))    # CVごとに，PGYの低い方と高い方のF，PGYを計算．ここでは80%PGYを計算しているが，0.8の数字を変えれば，他の値でも計算できる．

# 加入変動係数に対する平衡状態資源量の変化を見る図
pN <- ggplot(NY.curve)+geom_line(aes(x=F,y=Abundance,color=CV,size=I(3)))+theme_bw()+xlim(c(0,0.6))+theme(panel.grid.major = element_blank(),
 panel.grid.minor = element_blank(),
 axis.text=element_text(size=14),
 axis.title=element_text(size=18,face="bold"),
 legend.text=element_text(size=18),
 legend.title=element_text(face ="bold",size=18,hjust = 0)
 ) 

pN1 <- pN+stat_summary(data=MSY.stat,aes(x=Fmsy,y=Nmsy), fun.y=mean, colour = "red", size = 15, geom = "point", shape="*")    # MSYの点を加えてやっている

# 加入変動係数に対する平衡状態漁獲量の変化を見る図
pC <- ggplot(NY.curve)+geom_line(aes(x=F,y=Yield,color=CV,size=I(3)))+theme_bw()+xlim(c(0,0.6))+theme(panel.grid.major = element_blank(),
 panel.grid.minor = element_blank(),
 axis.text=element_text(size=14),
 axis.title=element_text(size=18,face="bold"),
 legend.text=element_text(size=18),
 legend.title=element_text(face ="bold",size=18,hjust = 0)
 )

pC1 <- pC+stat_summary(data=MSY.stat,aes(x=Fmsy,y=MSY), fun.y=mean, colour = "red", size = 15, geom = "point", shape="*")    # MSYの点を加えてやっている

# 資源量と漁獲量の図を並べる
grid.arrange(pN1,pC1,nrow=1)

# 加入変動の大きさによってMSYがどう異なるかが見られる．今回は，シンプルなモデルを使っているが，より複雑なモデルでも同じような結果になる．加入変動が大きくなると，目標資源量（SBmsy)が上がっていくことが見て取れる．これは，不確実性が大きくなると，それを考慮してより安全な目標を，という考え方の反映となっており，HSによる確率的MSYが持つ望ましい性質である（Beverton-HoltやRickerも同様の性質を持つが，HSはより明確である）．

## Transition in SBmsy
# 変動係数によるSBmsyの違いを図示する
# 変動係数が大きくなるとFが小さくなり目標値が上がる．これはprecautionary approachの考えと合致している

# 40年間のデータ
nY <- 40

rSB <- runif(nY,0,2000)    # 親魚量のデータを作ってやっている
z <- rnorm(nY)     # 加入変動の元になるものを作ってやる
SRrel <- NULL    # 結果をとっておく変数
k <- 1

# 加入変動係数を変えた場合のMSYを図示
for (sigma in c(0,0.2,0.4,0.6,0.8,1)){    # 加入変動（CV）の値の設定
  tR <- HS(rSB)     # 再生産関係から加入量を計算
  rR <- tR*exp(sigma*z)    # 上の加入量の誤差を加えている
  SRrel <- rbind(SRrel,cbind(CV=sigma,SB=rSB,R=rR,tR=tR,SBmsy=MSY.stat[k,4]*(1-MSY.stat[k,2])))     # 上で計算したMSYの結果から，MSYの親魚量をまとめたやっている
  k <- k + 1
}

SRrel <- as_tibble(SRrel)    # tibbleはtidyverseで使うデータ形式．dataframeよりグレードアップしている．dataframeでも良いが，ちょっと使ってみただけ．
colnames(SRrel) <- c("CV","SB","R","tR","SBmsy")

# 再生産データや曲線を描いた上で，そこにSBmsyの線をつけた．データの変動（加入量の変動）が大きくなると，SBmsyが大きくなっていくのが見られる．これは（当然ではあるが）決定論的なMSYにはない考え．単純な考え方ではあるが，持続的利用という観点で重要である．
ggplot(SRrel,aes(x=SB,y=R))+geom_point()+geom_line(aes(x=SB,y=tR))+geom_vline(aes(xintercept=Nmsy*(1-Fmsy),colour="red"),data=MSY.stat,linetype=2,size=1.5)+theme_bw()+theme(legend.position="none",text=element_text(size=16))+facet_wrap(~CV)

## Robustness of HS

# 100年間のデータを作って，データ数を増やして行った場合のHS再生産曲線やBH再生産曲線のbパラメータの変化を見てやる．HS再生産関係（頑健推定も利用）は，極端に大きなbパラメータとならず，かなり安定していることが見てとれる．本来は，SBmsyの変化を見たほうが良いかもしれないが，簡単なので再生産関係パラメータbを見てやっている．

source("SRfunctions.R")   # 再生産関係をフィットする関数を読み込む

nY <- 100    # シミュレーションの設定

set.seed(123)     # 乱数の設定．この値を設定しておけば，同じシミュレーション結果を再現できる．
rSB <- runif(nY,0,1500)    # 親魚量を作ってやる
z <- rnorm(nY)    # 加入変動の元となる変数
sigma <- 0.4     #  加入の変動係数=40%
xSB <- seq(0,1500,by=10)    # 再生産関係の予測値を描くためにSBのダミー変数を作っておく

tR <- BH(rSB)     # 真の再生産関係はBHなので，推定モデルとしてはBHが正しい．正しいモデルを使うことが常に良いわけではない，ということを見ようという意図．
rR <- tR*exp(sigma*z)     # 真の加入量に誤差を付す
SRdat <- as_tibble(cbind(CV=sigma,SB=rSB,R=rR,tR=tR))   #　再生産関係推定用のデータセットを作る

# SSBが350以下のデータの場合には，直線的な再生産関係のように見えるため，HS再生産曲線をフィットすると最大値が折点（b）となる．（正しい再生産関係である）BHを使うとbはすごく大きくなってしまう．HSを使って折点が最大値になる場合に対しておかしいという批判があるが，言われるほど不合理な考えではないのではないか？むしろ，なにがなんでも親魚量の観測範囲内に折点やSBmsyがあると考えないといけないという方に，どういう根拠があってそう考えるのか問いたい．この後のシミュレーションで，データが増えていくと，HSのb推定が正しい方に近づいていくのが見てとれる．
SRdat1 <- SRdat %>% filter(SB < 350)
HS1 <- estSR(list(R=SRdat1$R,SSB=SRdat1$SB),type="L1",SR="HS")
BH1 <- estSR(list(R=SRdat1$R,SSB=SRdat1$SB),type="L2",SR="BH")
p1 <- ggplot(SRdat1)+geom_point(aes(x=SB,y=R),size=1.5)+geom_vline(aes(xintercept=HS1$pars[2]),linetype=2,size=1.5,color="red")+geom_vline(aes(xintercept=BH1$pars[2]),linetype=2,size=1.5,color="blue")+theme_bw()+xlim(0,1600)+ylim(0,1200)+geom_line(aes(x=SB,y=R),data=as.data.frame(cbind(SB=xSB,R=pred.SR(HS1,xSB))),col="red",linetype=2,size=1.5)+geom_line(aes(x=SB,y=R),data=as.data.frame(cbind(SB=xSB,R=pred.SR(BH1,xSB))),col="blue",linetype=2,size=1.5)

# SSB < 800以下のデータに限定するとHSの折点の位置が少し右に動き，今度は観測範囲内に折点が出現するようになった．BHのbも見えるようになったが，観測範囲よりも大きく，真値よりもまだ大きい．
SRdat2 <- SRdat %>% filter(SB < 800)
HS2 <- estSR(list(R=SRdat2$R,SSB=SRdat2$SB),type="L1",SR="HS")
BH2 <- estSR(list(R=SRdat2$R,SSB=SRdat2$SB),type="L2",SR="BH")
p2 <- ggplot(SRdat2)+geom_point(aes(x=SB,y=R))+geom_vline(aes(xintercept=HS2$pars[2]),linetype=2,size=1.5,color="red")+geom_vline(aes(xintercept=BH2$pars[2]),linetype=2,size=1.5,color="blue")+theme_bw()+xlim(0,1600)+ylim(0,1200)+geom_line(aes(x=SB,y=R),data=as.data.frame(cbind(SB=xSB,R=pred.SR(HS2,xSB))),col="red",linetype=2,size=1.5)+geom_line(aes(x=SB,y=R),data=as.data.frame(cbind(SB=xSB,R=pred.SR(BH2,xSB))),col="blue",linetype=2,size=1.5)

# SB < 1200，SB < 1500と変えていくと，HSは比較的安定したbだが，BHのbは結構動く．この場合，真はBHであるので，HSが正しいわけではない．しかし，管理に対する目標値として使用するのに，どちらが良いかという話になるであろう．絶対的に正しいことをしようとしているわけではなく，不確実性に強い運用しやすい目標値の決め方は？という視点がある．このあたりは「資源評価」と「資源管理」の違いということもあるだろう．正しさを追求しすぎてしまうことは，時に「正しくない」....
SRdat3 <- SRdat %>% filter(SB < 1200)
HS3 <- estSR(list(R=SRdat3$R,SSB=SRdat3$SB),type="L1",SR="HS")
BH3 <- estSR(list(R=SRdat3$R,SSB=SRdat3$SB),type="L2",SR="BH")
p3 <- ggplot(SRdat3)+geom_point(aes(x=SB,y=R))+geom_vline(aes(xintercept=HS3$pars[2]),linetype=2,size=1.5,color="red")+geom_vline(aes(xintercept=BH3$pars[2]),linetype=2,size=1.5,color="blue")+theme_bw()+xlim(0,1600)+ylim(0,1200)+geom_line(aes(x=SB,y=R),data=as.data.frame(cbind(SB=xSB,R=pred.SR(HS3,xSB))),col="red",linetype=2,size=1.5)+geom_line(aes(x=SB,y=R),data=as.data.frame(cbind(SB=xSB,R=pred.SR(BH3,xSB))),col="blue",linetype=2,size=1.5)

SRdat4 <- SRdat %>% filter(SB < 1500)
HS4 <- estSR(list(R=SRdat4$R,SSB=SRdat4$SB),type="L1",SR="HS")
BH4 <- estSR(list(R=SRdat4$R,SSB=SRdat4$SB),type="L2",SR="BH")
p4 <- ggplot(SRdat4)+geom_point(aes(x=SB,y=R))+geom_vline(aes(xintercept=HS4$pars[2]),linetype=2,size=1.5,color="red")+geom_vline(aes(xintercept=BH4$pars[2]),linetype=2,size=1.5,color="blue")+theme_bw()+xlim(0,1600)+ylim(0,1200)+geom_line(aes(x=SB,y=R),data=as.data.frame(cbind(SB=xSB,R=pred.SR(HS4,xSB))),col="red",linetype=2,size=1.5)+geom_line(aes(x=SB,y=R),data=as.data.frame(cbind(SB=xSB,R=pred.SR(BH4,xSB))),col="blue",linetype=2,size=1.5)

# グラフを並べる
grid.arrange(p1,p2,p3,p4,nrow=2)

## What happens if environment controls all?
# もし環境だけで増加率が変わるような資源があるとしたら，条件を何も与えなければ，どこまでも増加していくか，絶滅するかしかなく，そもそもそんな資源に持続的利用を考える意味があるか？ということになる．少なくとも漁獲がない（F=0）ときに持続的であるためには，死亡と増加（加入）のバランスをとらないといけない．しかし，そうしても密度効果がないならば（ちょっとでも漁獲すれば絶滅に向かうため），持続的利用を考えることは難しい．だとしたら，密度効果があるとかないとかいう議論に意味があるだろうか....？そもそも我々は持続的利用が可能な資源（ではあるが不確実性が大きく真実が分かりにくい）に対する持続的利用は何か？という問題を考えているのではないだろうか？

# 環境影響が大きいので，それが魚の個体群動態を支配しているということがしばしば言われる．もちろん環境影響は大きい．しかし，たとえば今あるデータからMSYが推定できないということが，密度効果がないとか持続的漁獲（MSY）がない，ということとはイコールではない．持続的漁獲を考えるなら，持続的個体群を考えるべきであるというのは，一見してトートロジー的な話であり，机上の空論的に感じられるかもしれない．しかし，持続的でない個体群のシミュレーションモデルを作って，それの持続的管理を行う，という問題だとして，プログラムを作ってみれば，その無意味さに気づくことができるだろう．MSYがあるかないかなどは愚問であり，その意味では必ず「ある」のである（単純なMSYがあると言っているのではないことに注意）．もし，環境（だけ）が重要であるという人がいれば，あなたにとって持続的漁獲とは何か？と聞いてみて欲しい.... そもそも持続的でないもの（持続性を前提としないもの）を考えていて，それの持続的利用を考えているというなら，問題設定自体が間違っている，と考えられるでしょう．このことは，多くの人にとっては馬鹿馬鹿しい話であるかもしれないが，しかし，我が国では，そうしたことが割と受け入れられているのが現状である．環境の影響が大きいことは，資源管理不要となるのではなく，むしろ資源管理の必要性がずっと大きくなるということである（なぜなら環境の不確実性が大きいほど持続的利用は難しくなるから）．少なくとも，環境変動 vs. MSYというような話をしているのではなく，環境変動を取り込んだMSYとは？という話に移っているということを理解していただきたいところである．（ここでしたいことは，この考えを押し付けたいというような意図ではなく，何をどのように考えないといけないのかという理屈をこうした簡単な思考実験などを通して考えてみようというようなこと）

set.seed(1)
nS <- 50     # 増加率のパターンを50種類（1年続くレジーム，2年続くレジーム，...，50年続くレジーム）作る
a <- runif(nS,0.1,0.9)    #  増加率を0.1から0.9までの中で50個作る（上の年数のレジームに対応する増加率）

nY <- 1000    # 1000年間シミュレーションする

period <- rmultinom(1,nY,rep(1,nS))   #  1000年間における50の増加率パターンの出方の出現頻度を与えている
a.r <- a[unlist(lapply(1:nS, function(i) rep(i,period[i,])))]    # periodで与えたレジームが1000年間ランダムに出現する．
qplot(1:nY,a.r,geom="line",xlab="Year",ylab="a")+theme_bw()

N <- numeric(nY+1)    # シミュレーションした個体数を入れる
N[1] <- 1000    # 最初の個体数を1000にする
enN <- NULL
for (Fr in c(0.0,0.01,0.02)){  # ちょっとでも漁獲圧を大きくすると絶滅していってしまう
for (i in 1:nY){
S <- 1-a.r[i]     # Fr=0にすれば，持続的になる．この仮定がないと持続的にならず，無限に増加するか，絶滅するかになってしまう．
N[i+1] <- S*N[i]*(1-Fr)+a.r[i]*N[i]*(1-Fr)
}
enN <- rbind(enN, cbind(Fr,1:(nY+1),N))
}
enN <- data.frame(enN)
colnames(enN) <- c("Fr","Year","N")
enN[,1] <- as.factor(enN[,1])

# F=0, 0.01, 0.02による個体数の変化をプロット
ggplot(enN)+geom_line(aes(x=Year,y=N,color=Fr),data=enN,size=3)+theme_bw()+theme(panel.grid.major = element_blank(),
 panel.grid.minor = element_blank(),
 axis.text=element_text(size=14),
 axis.title=element_text(size=18,face="bold"),
 legend.text=element_text(size=18),
 legend.title=element_text(face ="bold",size=18,hjust = 0)
 ) 
 
# ここでのメッセージは，世の中のレジームシフト信仰のようなものがどういう意味を持っているか考えていただきたいというようなこと．レジームシフトがない，と言っているわけではない．そのような変動環境のもとで持続的であるように適応した個体群を考えないと意味がないのではないですか？ということ．持続性のない資源もあるかもしれないが（今，まさに絶滅に向かっている），少なくとも前提としてそういうものを（持続的に）管理しようとしているわけではない．「MSY理論は嘘だ」というような話があるが，それを鵜呑みにする前に，少し考えてみていただきたい．我々ももちろん特定の（古典的な）MSYのような量があると言っているわけではない．しかし，MSYは有用な道具（のひとつ）であるという点では，Daniel PaulyやAndre Puntと同意見である．我が国でよく引用されるレジームシフトを重視してMSYを批判するような言説にどういった意味があるのだろうか？レジームシフトがある中で，実用的なMSYはどういうものだろうか？という発想になぜいかないのだろうか....？（MSY批判にやっきになるだけではなく，きちんとした対案を出して，それをMSEなどで比較して良いHCRを選ぶことが重要なのでは....？）

## Is a complex model better than a simple model?
# HS再生産曲線は再生産データにうまくあてはまっていないではないかという批判がある．しかし，再生産曲線は非常にノイズが大きいのが普通で，それにあてはまる曲線を見つけるのは容易ではない．そのような大きなノイズの下でも頑健に管理できる方法は何か？と考えるべきであろう．複雑なモデルが良いわけではない．欲しいのは，予測力の高いモデルである（参考：資源管理研修会上級編の機械学習のスライド）

set.seed(13)
nY <- 50    # 50年間のデータ
rSB <- runif(nY,0,1500)    # 親魚量を作る
z <- rnorm(nY)    # 加入尾数の誤差のベース
sigma <- 0.4    # 加入変動の大きさ

tR <- BH(rSB)    # 真の再生産関係はBH
rR <- tR*exp(sigma*z)    # 加入変動を付与
SRdat <- as_tibble(cbind(CV=sigma,SB=rSB,R=rR,tR=tR))    # 再生産関係推定のためのデータセット

# SB < 1300をモデルフィットに使い，SB >= 1300を予測する問題を考える
train.SR <- SRdat %>% filter(SB < 1300)   # モデルフイットのためのデータ
test.SR <- SRdat %>% filter(SB >= 1300)   # 親魚量が1300以上のときをテストデータにする（外挿してうまく推定できるかみたいので）

cmod <- smooth.spline(train.SR$SB,log(train.SR$R),spar=0.1)    # スプライン曲線で複雑なのをあてはめてみた（sparで複雑さを調整している）

fit.HS <- estSR(list(R=train.SR$R,SSB=train.SR$SB),type="L1",SR="HS")   # SRfunctionsに入っている再生産関係のあてはめ関数．type="L1"で頑健推定（絶対値最小化），SR="HS"でHockey-Stickを指定している

# 図を見ると，手元にあるデータにあまりあてはまっていないようなHSの方が，うまく予測できているように見える．加入の変化を忠実に追うようなモデルが良いわけではなく，予測困難な加入の変化に頑健なモデルを選ぶことはより重要である．

plot(train.SR$SB, train.SR$R, xlim=c(0,1500), ylim=c(0,1500), cex=2, pch=16, xlab="SB", ylab="R")
lines(predict(cmod,sort(train.SR$SB))$x, exp(predict(cmod,sort(train.SR$SB))$y), col="orange", lwd=3)
lines(sort(train.SR$SB), pred.SR(fit.HS,sort(train.SR$SB)), col="red", lwd=2)

points(test.SR$SB, test.SR$R, col="blue", cex=2, pch=16)

lines(predict(cmod,sort(test.SR$SB))$x, exp(predict(cmod,sort(test.SR$SB))$y), col="orange", lwd=4)

lines(sort(test.SR$SB), pred.SR(fit.HS,sort(test.SR$SB)), col="red", lwd=4)

# HS再生産曲線がデータにあてはまっていない．だから有用ではない，信頼できない，という批判がある．そして加入データによく適合するモデルが良いという話もあるが，必ずしもそうではない．我々がしたいことは，未知のデータの予測（将来予測）である．したがって，手元のデータによく適合するモデルが良いモデルというのを妄信してはならない．手元のデータにあまりに当てはまりすぎるものはむしろ信用ならないのである．このことは，統計学でバイアス・分散トレードオフと呼ばれている．加入変動を環境情報で予測しようという研究もあるが，それが非常に成功したという例はいまだ少ない．特に長期的に予測することはかなり難しくほぼ不可能であろう（短期的な予測には希望がないわけではない）．そのような前提のもとで，長期的に安定的な漁業（つまり，持続的な漁業）をどう実現しようか，と考えるべきではないだろうか？

## MSY under regime shift
# レジームシフトのもとで新ルール（のMSY）はどのような振る舞いをするだろうか？レジームシフトがあるからといって，全然駄目にはならない．なぜなら，不確実性に強いMSYの考え（不確実性が大きければ保守的になる）とPretty Good Yield（多少保守的な（小さい）漁獲圧になっても漁獲量としてはそんなに損しない）の考えがあるから．

# estMSY function
# HS再生産曲線をフィットした結果を使って資源量や漁獲量をシミュレーションする関数

estMSY <- function(Fr,Dat,seed=1){
  set.seed(seed)
  resSR <- estSR(list(SSB=Dat$SB,R=Dat$R),type="L1",SR="HS")    # 再生産関係をfitする
  R <- N <- numeric(51)
  N[1] <- resSR$pars[2]      # bをNの初期値とする
  R[1] <- N[1]*(1-S)     # N=S*N+Rから
  z <- rnorm(50)      # 加入変動に関係する誤差の変動項
    
  for (i in 1:50){
    R[i+1] <- HS(N[i]*(1-Fr),a=resSR$pars[1],b=resSR$pars[2])*exp(resSR$pars[3]*z[i]-resSR$pars[3]^2/2)    # 再生産関係モデルから推定されたパラメータで加入量を計算
    N[i+1] <- S*N[i]*(1-Fr)+R[i+1]     # 個体群動態モデル
  }
  
  list(N=N, SB=N*(1-Fr), C=N*Fr)     # 資源量，親魚量，漁獲量
}

# regime shift simulation
# レジームシフト資源を作りだす

# a1[1], b1[1]は高レジームのパラメータ，a1[2], b1[2]は低レジームのパラメータ．R0が10倍違う
a1 <- c(1.4,0.7)
b1 <- c(1000,200)

# 100年間シミュレーションする
nY <- 100

set.seed(1)

# レジームシフトのパターン（高→低）は100年で4回繰り返される．高，低の出現確率は半々
rs <- rbinom(4,25,0.5)     #  25年の高低1セットを4回繰り返すので，25年のうち高レジームの期間を乱数で発生する

# レジームシフトの状態を記録する
rs.state <- unlist(map(1:4, function(i) c(rep(1,rs[i]),rep(2,len=25-rs[i]))))    # 25年間で高レジームと低レジームがおよそ半分ずつで交替する

# 100年間シミュレーションする
N <- R <- numeric(nY+1)
N[1] <- a1[1]*b1[1]/(1-S)    # 高レジームのときのN0からスタート
R[1] <- a1[1]*b1[1]      # 高レジームのときのR0
z <- rnorm(nY)
sigma <- 0.2    # それぞれのレジームにおける加入変動
Fr <- rlnorm(nY,log(0.6*(1-1/(a1[rs.state]+S))),0.2)     # 漁獲率は決定論的Fmsyの60%の平均値で確率的にふらせる．どういうFにするか，ここでは適当に決めている．

for (i in 1:nY){
  R[i+1] <- HS(N[i]*(1-Fr[i]),a=a1[rs.state[i]],b=b1[rs.state[i]])*exp(sigma*z[i]-sigma^2/2)    # レジームシフトによって変動する加入
  N[i+1] <- S*N[i]*(1-Fr[i])+R[i+1]    # 個体群動態
}

# シミュレーションした個体群の様子
qplot(x=1:100, y=N[1:100], size=I(5), color=factor(rs.state))+guides(color=FALSE)+theme(axis.text=element_text(size=14), axis.title=element_text(size=18,face="bold"))+theme_bw()+xlab("Year")+ylab("Abundance")

# シミュレーションした個体群をレジームシフトの状態で色分けしている．ここで注意すべきことは，個体群は時系列でつながっているので，たとえ再生産パラメータが変化したとしても，その個体群への影響は遅れて現れるということである（ということは，レジームの判別はより難しいことになる）．

# 高レジーム，低レジームそれぞれのMSYを計算してやる
Res.rsMSY <- list()
rsNY.curve <- NULL
k <- 1
Fr.range <- seq(0,0.99,by=0.01)
for (k in 1:2){
  Res.rsMSY[[k]] <- map(Fr.range, Sim.msy, a=a1[k], b=b1[k], sigma=sigma)

  rsNY.curve <- rbind(rsNY.curve, cbind(k,Fr.range,map_dbl(1:100,function(i) mean(Res.rsMSY[[k]][[i]]$N)),map_dbl(1:100,function(i) mean(Res.rsMSY[[k]][[i]]$C))))
}

rsNY.curve <- as.data.frame(rsNY.curve)
colnames(rsNY.curve) <- c("RS","F","Abundance","Yield")
rsNY.curve[,1] <- factor(rsNY.curve[,1])

rsMSY.stat <- rsNY.curve %>%
  group_by(RS) %>%
  summarise(Fmsy=F[which.max(Yield)],MSY=max(Yield),Nmsy=Abundance[which.max(Yield)])

# Fや再生産パラメータが毎年変わってもシミュレーションできる関数
rsSim.msy <- function(Fr=0, Y=100, sigma=0.2, a=1.2, b=1000, seed=1){
  set.seed(seed)
  N <- C <- numeric(Y)
  N[1] <- Nmsy
  z <- rnorm(Y)
  
  for (i in 1:Y){
    C[i] <- N[i]*Fr[i]
    N[i+1] <- pmax(S*(N[i]-C[i])+SRfunc(N[i]-C[i],a=a[i],b=b[i],SR=SR)*exp(sigma*z[i]-sigma^2/2),0)
  }
  
  list(C=C,N=N)
}

# 再生産曲線推定のためのデータを作る
Dat <- data.frame(R=R[-1], SB=N[1:nY]*(1-Fr))

# レジームシフトの出現パターンは正確に知っており，それに基づいて“正しい”Fmsyを使う．Fmsyはデータから推定する
Fr.range1 <- seq(0,0.9,by=0.1)

Fmsy1 <- Fr.range1[which.max(map_dbl(Fr.range1, function(Fr) mean(estMSY(Fr,Dat[rs.state==1,])$C)))]   # MSYは再生産関係のデータから推定するが，レジームのパターンごとに推定してやる．ここでは，高レジームの場合．大雑把に推定して，最適値の周辺を判断し，最適化関数を使って最適値を求める．
res.msy1 <- optimize(function(x) -mean(estMSY(x,Dat[rs.state==1,])$C),c(Fmsy1-0.1,Fmsy1+0.1))
Fmsy1 <- res.msy1$minimum
MSY1 <- res.msy1$objective*(-1)

Fmsy2 <- Fr.range1[which.max(map_dbl(Fr.range1, function(Fr) mean(estMSY(Fr,Dat[rs.state==2,])$C)))]   # 低レジームの場合のMSY
res.msy2 <- optimize(function(x) -mean(estMSY(x,Dat[rs.state==2,])$C),c(Fmsy2-0.1,Fmsy2+0.1))
Fmsy2 <- res.msy2$minimum
MSY2 <- res.msy2$objective*(-1)

rsFmsy <- c(Fmsy1,Fmsy2)
rsFmsy <- rsFmsy[rs.state]

Res.rsMSY <- rsSim.msy(rsFmsy, a=a1[rs.state], b=b1[rs.state], sigma=sigma)     # レジーム状態に正確に対応したFmsyで漁獲した場合の資源量や漁獲量を計算．

# レジームシフトモデルによるNmsyとMSYの平均値，中央値
rsN.mean <- mean(Res.rsMSY$N)
rsC.mean <- mean(Res.rsMSY$C)
rsN.median <- median(Res.rsMSY$N)
rsC.median <- median(Res.rsMSY$C)

# 新ルールに従って全体でひとつのFmsyを推定してやる
Fmsy <- Fr.range1[which.max(map_dbl(Fr.range1, function(Fr) mean(estMSY(Fr,Dat)$C)))]      # 再生産関係データ全体を使ってMSYを推定する．最初は大雑把に推定し，optimizeを使って細かい推定をする．
res.msy <- optimize(function(x) -mean(estMSY(x,Dat)$C),c(Fmsy-0.1,Fmsy+0.1))
Fmsy <- res.msy$minimum
MSY <- res.msy$objective*(-1)
Res.newMSY <- rsSim.msy(rep(Fmsy,nY), a=a1[rs.state], b=b1[rs.state], sigma=sigma)      #　全体で推定したFmsyでレジームシフトで変わる資源を漁獲していった結果

# 新ルールによるNmsy，MSYの平均値，中央値
newN.mean <- mean(Res.newMSY$N)
newC.mean <- mean(Res.newMSY$C)
newN.median <- median(Res.newMSY$N)
newC.median <- median(Res.newMSY$C)

# レジームシフトと新ルールのMSYを比較
Res.RS <- data.frame(Type=c("C:mean","C:median","N:mean","N:median"),Ratio=c(newC.mean/rsC.mean,newC.median/rsC.median,newN.mean/rsN.mean,newN.median/rsN.median))

ggplot(Res.RS)+geom_bar(aes(x=Type,y=Ratio),fill=3,stat="identity")+geom_hline(aes(yintercept=1),color="red",size=3,linetype=2)+coord_flip()+theme_bw()+theme(axis.text=element_text(size=18),
 axis.title=element_text(size=18,face="bold"),
 legend.text=element_text(size=18),
 legend.title=element_text(face ="bold",size=18,hjust = 0))

# レジームシフトによる変化が大きいので，それを考慮しない管理はうまくいかないような印象を持つが，確率的MSYにより変動が大きい場合は（レジーム変動も含めて誤差とするので大きな変動になる）保守的な目標値になるので，持続性が損なわれることはない．しかし，漁獲率が低くなることにより，漁獲の効率が悪くなると想像されるが，結構保守的にしてもPretty Good Yieldによりそこまで漁獲は減らない．結局，レジームを考慮しないでも，そこまで非効率な漁獲にならず，持続的な漁獲を行えることになる．しかし，だから，レジームシフト資源に対して，新ルールでも絶対に大丈夫ということではない．レジームシフトの程度によってきめ細かい対応が必要になるだろう．ここでは，我々が単純なMSYで想像するほど単純に駄目ではなく，新しいMSYの考えを活用することにより，レジームシフト資源にもうまく対応できるだろう，ということが言いたいことである．ここでもレジームシフトのような状態にあるというのを現在のデータを使って示すことではなく，レジームシフトを精度良く将来予測ができるかどうかが重要である．将来予測が高精度で行えるならそれを（より効率の良い持続的漁業のために）活用する余地は十分にある．しかし，将来予測が精度良くできないなら，それは不確実性としてとらえ，そのような不確実性に頑健な漁獲方式を採用すべきである．

## Pretty Good Yield
# PGYのときの資源量がどんな感じか計算して図示
# PGYの考え方はアメリカの資源管理や我が国の新しい資源管理に大きな影響を与えている

NY.curve1 <- NY.curve %>% group_by(CV) %>% mutate(isPGY = Yield >= 0.8*max(Yield), relAbund=Abundance/max(Abundance))

p1 <- ggplot(NY.curve1)+geom_point(aes(x=F, y=relAbund, color=isPGY), size=2)+facet_wrap(~CV)+theme_bw()+xlim(0,0.6)
p1

## Abundance at MSY
# Fmsyで漁獲し続けた場合の資源量が確率変動のもとでどれだけ変わるかを図示
# Fmsyは平均的に漁獲を最大化するものなので，その変動が大きい場合，厳しいリスクを避けられない場合があるかもしれない

plot(Res.sMSY[[5]][[30]]$N,ylim=c(0,5000),pch=16,xlab="Year",ylab="Abundance",col="blue",cex=1.2,main="CV = 0.8")
abline(h=MSY.stat[5,4],col="red",lwd=2,lty=2)

# MSYは一定のFで平均漁獲量が最大になるものとしているので，そのFで獲り続けたとしても，それがいつでもベストではない．加入変動以外の不確実性もあることを考えれば，Fmsyに少し手を加えて，より良いルールのもとに漁獲を行うのが良いだろうと想像できる．

## HCR
# Harvest Control Ruleの簡単なMSE的なシミュレーション

# 実際にはもう少し複雑なことをしているが，大体こんなようなシミュレーションで，HCRの良さを評価している．このようなシミュレーションによって頑健性を確認した上で，HCRを決めるのが基本的なことではある．

# 新ルールのHCR
HCR <- function(x,beta=0.8,Fmsy=0.2,SBlim=0,SBban=0) beta*Fmsy*pmax(pmin((x-SBban)/(SBlim-SBban),1),0)   # 新ルールのHCR（の簡易版）

# pre-management data
# 管理開始前のデータを作成してやる．開始前に40年のデータがあるとする

  set.seed(1)
  N <- R <- numeric(40)
  N[1] <- N0
  R[1] <- R0
  preY <- 40-1
  postY <- 50    # 将来予測は50年
  Fr <- rlnorm(preY,log(as.numeric(PGY.stat[1,4])),0.2)    # 例によって，歴史的な漁獲率を適当に与えている．40年後の資源量が20%N0になるように調整したりもできるが，ここではコードを簡単にするため適当な値を入れている．
  z <- rnorm(preY+postY)    # 加入変動の変数は管理前後の年数分だけ作っておく
  sigma <- 0.6     # 加入変動の大きさを60%に設定
  
  trueSBmsy <- as.numeric(MSY.stat[4,4]*(1-MSY.stat[4,2]))       # 真のSBmsy
  trueFmsy <- as.numeric(MSY.stat[4,2])    # 真のFmsy
  
  for (i in 1:preY){
    R[i+1] <- SRfunc(N[i]*(1-Fr[i]),SR=SR)*exp(sigma*z[i])     # 加入量
    N[i+1] <- S*N[i]*(1-Fr[i])+R[i+1]    # 個体群動態
  }

  preN <- N[1:preY]
  preDat <- data.frame(N=preN,SB=preN*(1-Fr),C=preN*Fr,R=R[2:(preY+1)])     # 管理開始前のデータを作ってやる

  postDat <- preDat

# forecasting
# 将来予測をする

Fr.range <- seq(0,0.9,by=0.1)
beta <- 0.8     # HCRの中の調整係数

for (k in 1:postY){
  # 資源量を更新
  pre.N <- last(postDat$N)    # データセットの最終年の資源量
  pre.SB <- last(postDat$SB)    # データセットの最終年の親魚量
  
  post.R <- max(SRfunc(pre.SB,SR=SR)*exp(sigma*z[preY+k]-sigma^2/2),0)     # 親魚量から加入量を計算
  post.N <- max(S*pre.SB+post.R,0)    # 資源量を更新
 
  if (k %% 5 == 1){   # 5年に1回管理基準値を更新
    Fmsy <- Fr.range[which.max(map_dbl(Fr.range, function(Fr) mean(tail(estMSY(Fr,postDat)$C,10))))]     # シミュレーションの最後10年間の平均漁獲量が最大になるものをMSYと定義している（結果を安定させるため）．まず0.1刻みぐらいでおおざっぱに計算する
    res.msy <- optimize(function(x) -mean(tail(estMSY(x,postDat)$C,10)),c(Fmsy-0.1,Fmsy+0.1))     # 細かい数値まで最適化してFmsyを計算
    Fmsy <- res.msy$minimum
    MSY <- res.msy$objective*(-1)
    Flim <- uniroot(function(x) 0.6*MSY-mean(tail(estMSY(x,postDat)$C,10)),c(Fmsy,0.99))$root     # Flimitを計算
    Fban <- uniroot(function(x) 0.1*MSY-mean(tail(estMSY(x,postDat)$C,10)),c(Fmsy,0.99))$root     # Fbanを計算
    SBmsy <- mean(tail(estMSY(Fmsy,postDat)$SB,10))   # SBmsyを計算
    SBlim <- mean(tail(estMSY(Flim,postDat)$SB,10))   # SBlimを計算
    SBban <- mean(tail(estMSY(Fban,postDat)$SB,10))   # SBbanを計算
  }
  
  ABC <- HCR(pre.SB,beta=beta,Fmsy=Fmsy,SBlim=SBlim,SBban=SBban)*post.N    # HCRに基づきABCを計算
  
  post.SB <- max(post.N-ABC,0)     #  親魚量を更新
  post.C <- ABC    # 漁獲量を更新
  
  postDat <- rbind(postDat, c(post.N, post.SB, post.C, post.R))
}

# 真のSBmsy，Fmsyに対して実際のSB，Fがどうかというプロットを作ってやってHCRによる管理の性能を評価する
performance.stat <- data.frame(management.period=1:(preY+postY),SBoverSBmsy=postDat$SB/trueSBmsy, FoverFmsy=(postDat$C/postDat$N)/trueFmsy)

ggplot(performance.stat)+geom_point(aes(x=SBoverSBmsy, y=FoverFmsy, color=management.period), size=3)+theme_bw()
