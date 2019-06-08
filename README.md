# Riemann solver
Riemann問題の厳密解を求めるコードです。

## Riemann問題
![shocktube.png](https://raw.githubusercontent.com/mitchiinaga/riemann_solver/master/images/riemann_image.png)

時刻t=0で、x=0を境に右側にU_R、左側にU_Lの物理量を持った1次元非粘性圧縮性流体が存在しているとします。
この流体の物理量の時間発展は、非線形方程式を解くことによって、厳密に求めることができます。
これがRiemann問題と呼ばれている初期値問題です。

詳細はToro (2009)などを参照してください。

## パラメータ
jsonファイルを引数に指定してパラメータを設定します。パラメータファイルは以下のように記述します。

```json:shocktube.json
{
    "left"  : {
        "velocity" : 0.0,
        "density"  : 1.0,
        "pressure" : 1.0
    },
    "right" : {
        "velocity" : 0.0,
        "density"  : 0.125,
        "pressure" : 0.1
    },
    "gamma"    : 1.4,
    "number"   : 100,
    "time"     : 0.2,
    "fileName" : "shocktube.dat"
}
```

それぞれのパラメータの意味は以下の通りです。
* left: x &lt; 0 の物理量。指定しなければならない物理量は速度(velocity)、密度(density)、圧力(pressure)の３つです。
* right: x &gt; 0 の物理量。
* gamma: 比熱比。
* number: 出力される点の数。
* time: 出力時刻。
* fileName: 出力ファイルの名前。

サンプルとして上記のshocktube.jsonに加え、123.json、blustwave.jsonの3種類を用意しています。
* shocktube.json: Sod Shock Tubeと呼ばれる、数値流体コードのテストによく使われる問題です(Sod 1978)。右側に衝撃波、左側に膨張波が現れます。
* 123.json: 123 problemやSj&ouml;green testなどと呼ばれる問題です(Einfeldt et al. 1991)。両側に膨張波が生じ、その内側が非常に低密度な領域となります。
* blustwave.json: 右側と左側の圧力比が100,000倍あり、非常に強い衝撃波が生じます。しかし、ランキン・ユゴニオの関係式から、密度の上昇は高々6倍程度 (比熱比が7/5の場合) となります。

## 計算例
shocktube.jsonを使用してSod Shock Tube問題を解いた結果を下に示します。
![shocktube.png](https://raw.githubusercontent.com/mitchiinaga/riemann_solver/master/images/shocktube_.png)

## 参考文献
* Toro, E. F. (2009). Riemann Solvers and Numerical Methods for Fluid Dynamics: A Practical Introduction (Vol. 53). Berlin, Heidelberg: Springer Berlin Heidelberg. https://doi.org/10.1007/b79761
* Sod, G. A. (1978). A survey of several finite difference methods for systems of nonlinear hyperbolic conservation laws. Journal of computational physics, 27(1), 1-31.
* Einfeldt, B., Munz, C. D., Roe, P. L., & Sj&ouml;green, B. (1991). On Godunov-type methods near low densities. Journal of Computational Physics, 92(2), 273–295. https://doi.org/10.1016/0021-9991(91)90211-3
