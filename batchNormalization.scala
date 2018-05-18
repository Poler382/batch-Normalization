import breeze.linalg._
import math._

object BN{
  def sub(y:Array[DenseVector[Double]],a:Array[Array[Double]] )={
    var d = Array.ofDim[Double](a.size,a(0).size)
    for(i <- 0 until a.size){
      for(j <- 0 until a(0).size){
        d(i)(j) = y(i)(j) - a(i)(j)
      }
    }
    d
  }

  def test()={
    val bn   = new BN(2,3)
    val test = Array(Array(1d,2d),Array(3d,4d),Array(5d,6d)) 

    val l = bn.forward(test)

    bn.backward(sub(l,test))
  }


}

class BN(val D:Int,val n: Int){
  var gamma = DenseVector.zeros[Double](n).map(a => 1d)
  var beta  = DenseVector.zeros[Double](n)
  var mu    = DenseVector.zeros[Double](D)
  var eps   = 0.00000001
  var sigma = DenseVector.zeros[Double](n)


  def forward(x:Array[Array[Double]]):Array[DenseVector[Double]] ={
    var out = new Array[DenseVector[Double]](x.size)
    
    //ave
    var sum =new Array[Double](D)
    for(i <- 0 until n){
      for(j <- 0 until D){
        sum(j) += x(i)(j)
      }
    }
    for(i <- 0 until D){
      mu(i) = sum(i) / x.size.toDouble 
    }    
    
    var xmu = Array.ofDim[Double](n,D)

    for(i <- 0 until n){
      for(j <- 0 until D){
        xmu(i)(j)= x(i)(j) - mu(j)
      }
    }

    var sigma = DenseVector.zeros[Double](D)
    sum =new Array[Double](D)
    for(i <- 0 until n){
      for(j <- 0 until D){
        sum(j) += xmu(i)(j) * xmu(i)(j)  
      }
    }
   
    for(i <- 0 until D){
      sigma(i) =1 / sqrt( sum(i) / n.toDouble)
    }

    for(i <- 0 until n){
      var xhat = DenseVector.zeros[Double](D)
      for(j <- 0 until D){
        xhat(j) = sigma(j) * xmu(i)(j) * gamma(i)+beta(i)
      }
      out(i) = xhat
    }
    out
  }

  def sumMatrix(in:Array[Array[Double]])={
    var m  = DenseVector.zeros[Double](D)
    var sum = new Array[Double](D)
    for(i <- 0 until n){
      for(j <- 0 until D){
        sum(j) += in(i)(j)
      }
    }
    for(i <- 0 until D){
      m(i) = sum(i) / in.size.toDouble 
    }

    m.toArray
  }
  def backward(d:Array[Array[Double]]) = {
    var d_beta  = sumMatrix(d)
    var d_gamma =



    d
  }

  def update()={

  }

  def reset()={



  }

}
