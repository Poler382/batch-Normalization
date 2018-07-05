class GroupNormalization (N:Int,C:Int,gr:Int){
  val rand = new util.Random(0)

  def forward(xs: Array[Array[Double]])={
    val chnsize = xs(0).size /C
    val bs = xs.size
    var myulist = new Array[Double](10)
    val x_mu = Array.ofDim[Double](xs.size,xs(0).size)
    var grlist = new Array[Array[Double]](gr)
    val sum = new Array [Double](gr)
    var ave2 = new Array [Double](gr)

    for (j <- 0 until xs(0).size ){
      for(i <- 0 until xs.size){
        sum(j/chnsize)+=xs(i)(j)
      }
    }

    for (j <- 0 until xs(0).size ){
      for(i <- 0 until xs.size){
        x_mu(i)(j) = xs(i)(j) - sum(j/chnsize) /(chnsize * bs)
      }
    }

    var x2 =x_mu.map(_.map(a => (a * a)/(chnsize * bs)))

    for (j <- 0 until xs(0).size ){
      for(i <- 0 until xs.size){
        ave2(j/chnsize)+=x2(i)(j) / (chnsize * bs)
      }
    }

    var sigma = ave2.map(a => 1d / math.sqrt(a+rand.nextGaussian()))

    println(sigma(0))
    for(i <- 0 until xs.size){
      for(j <- 0 until xs(0).size){
        println((i,j)+" : "+xs(i)(j)
                     +" -> "+x_mu(i)(j)
                     +" -> "+x2(i)(j))+" -> "+sigma(j/chnsize)

      }
      println()
    }


  }

  
  def cal_myu(x:Array[Double])=x.sum/x.size
  def cal_sigma(x:Array[Double],myu:Double)={
    val M = new Array[Double](x.size).map(_ => myu)
    val t = (x.zip(M).map{case (a,b) => a*a-2*a*b+b*b)}.sum )/ x.size
    
    math.sqrt(rand.nextGausssian+t)
  }
  def forward2(xs:Array[Array[Double]])={
    val chnsize = xs(0).size /C
    val bs = xs.size
    var myulist = new Array[Double](10)
    val x_mu = Array.ofDim[Double](xs.size,xs(0).size)
    var grlist = new Array[Array[Double]](gr)
    val sum = new Array [Double](gr)
    var ave2 = new Array [Double](gr)
    var channel_x = Array.ofDim[Double](gr,chnsize)

    for (j <- 0 until xs(0).size ){
      for(i <- 0 until xs.size){
        channel_x(j/chnsize)=xs(i)(j)
      }
    }


  }


}

object t{
  def main(args: Array[String]): Unit ={
    val xs = Array(
      Array(1d,2d,3d, 4d, 5d, 6d, 7d, 8d),
      Array(4d,5d,6d, 7d, 8d, 9d,10d,11d),
      Array(7d,8d,9d,10d,11d,12d,13d,14d)
    )
    val gn = new GroupNormalization(3,2,2)
    gn.forward(xs)

  }

}
