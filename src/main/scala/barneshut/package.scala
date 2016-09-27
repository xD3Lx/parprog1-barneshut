import common._
import barneshut.conctrees._

package object barneshut {

  class Boundaries {
    var minX = Float.MaxValue

    var minY = Float.MaxValue

    var maxX = Float.MinValue

    var maxY = Float.MinValue

    def width = maxX - minX

    def height = maxY - minY

    def size = math.max(width, height)

    def centerX = minX + width / 2

    def centerY = minY + height / 2

    override def toString = s"Boundaries($minX, $minY, $maxX, $maxY)"
  }

  /*
  mass = m_B + m_C + m_D + m_E
massX = (m_B * x_B + m_C * x_C + m_D * x_D + m_E * x_E) / mass
massY = (m_B * y_B + m_C * y_C + m_D * y_D + m_E * y_E) / mass
   */

  sealed abstract class Quad {
    def massX: Float

    def massY: Float

    def mass: Float

    def centerX: Float

    def centerY: Float

    def size: Float

    def total: Int

    def insert(b: Body): Quad
  }

  case class Empty(centerX: Float, centerY: Float, size: Float) extends Quad {
    def massX: Float = centerX
    def massY: Float = centerY
    def mass: Float = 0
    def total: Int = 0
    def insert(b: Body): Quad = Leaf(centerX, centerY, size, Seq(b))
  }

  case class Fork(
    nw: Quad, ne: Quad, sw: Quad, se: Quad
  ) extends Quad {
    val centerX: Float = nw.centerX + nw.size / 2
    val centerY: Float = nw.centerY + nw.size / 2
    val size: Float = nw.size + sw.size
    val mass: Float = nw.mass + ne.mass + sw.mass + se.mass
    val total: Int = nw.total + ne.total + sw.total + se.total
    val massX: Float = if (total == 0) centerX else getMassX(nw) + getMassX(ne) + getMassX(sw) + getMassX(se)
    val massY: Float = if (total == 0) centerY else getMassY(nw) + getMassY(ne) + getMassY(sw) + getMassY(se)

    def getMassX(q: Quad) = q.massX * q.mass / mass
    def getMassY(q: Quad) = q.massY * q.mass / mass

    def insert(b: Body): Fork = {
      val bx = b.x
      val by = b.y

      if ((bx <= centerX) && (by <= centerY)) {
        copy(nw = nw.insert(b))
      } else if ((bx <= centerX) && (by > centerY)) {
        copy(sw = sw.insert(b))
      } else if ((bx > centerX) && (by <= centerY)) {
        copy(ne = ne.insert(b))
      } else {
        copy(se = se.insert(b))
      }
    }
  }

  case class Leaf(centerX: Float, centerY: Float, size: Float, bodies: Seq[Body])
  extends Quad {
    val mass = tMass(bodies)
    val (massX, massY) = cMass(bodies, mass)
    val total: Int = bodies.size
    def insert(b: Body): Quad = {
      if (size > minimumSize) {
        val quadSize = size / 4

        val emptyFork = Fork(
          Empty(centerX - quadSize, centerY - quadSize, quadSize * 2),
          Empty(centerX + quadSize, centerY - quadSize, quadSize * 2),
          Empty(centerX - quadSize, centerY + quadSize, quadSize * 2),
          Empty(centerX + quadSize, centerY + quadSize, quadSize * 2)
        )

        (bodies :+ b).foldLeft(emptyFork) (_.insert(_))
      } else Leaf(centerX, centerY, size, bodies :+ b)
    }
  }

  private def tMass(bodies: Seq[Body]): Float = bodies.map(_.mass).sum
  private def cMass(bodies: Seq[Body], tMass: Float): (Float, Float) = {
    val sum = bodies
      .map { b => (b.x * b.mass, b.y * b.mass) }
      .fold((0f, 0f)) { case (acc, elem) => (acc._1 + elem._1, acc._2 + elem._2) }

    (sum._1 / tMass, sum._2 / tMass)
  }

  def minimumSize = 0.00001f

  def gee: Float = 100.0f

  def delta: Float = 0.01f

  def theta = 0.5f

  def eliminationThreshold = 0.5f

  def force(m1: Float, m2: Float, dist: Float): Float = gee * m1 * m2 / (dist * dist)

  def distance(x0: Float, y0: Float, x1: Float, y1: Float): Float = {
    math.sqrt((x1 - x0) * (x1 - x0) + (y1 - y0) * (y1 - y0)).toFloat
  }

  class Body(val mass: Float, val x: Float, val y: Float, val xspeed: Float, val yspeed: Float) {

    def updated(quad: Quad): Body = {
      var netforcex = 0.0f
      var netforcey = 0.0f

      def addForce(thatMass: Float, thatMassX: Float, thatMassY: Float): Unit = {
        val dist = distance(thatMassX, thatMassY, x, y)
        /* If the distance is smaller than 1f, we enter the realm of close
         * body interactions. Since we do not model them in this simplistic
         * implementation, bodies at extreme proximities get a huge acceleration,
         * and are catapulted from each other's gravitational pull at extreme
         * velocities (something like this:
         * http://en.wikipedia.org/wiki/Interplanetary_spaceflight#Gravitational_slingshot).
         * To decrease the effect of this gravitational slingshot, as a very
         * simple approximation, we ignore gravity at extreme proximities.
         */
        if (dist > 1f) {
          val dforce = force(mass, thatMass, dist)
          val xn = (thatMassX - x) / dist
          val yn = (thatMassY - y) / dist
          val dforcex = dforce * xn
          val dforcey = dforce * yn
          netforcex += dforcex
          netforcey += dforcey
        }
      }

      def traverse(quad: Quad): Unit = (quad: Quad) match {
        case Empty(_, _, _) =>
          // no force
        case Leaf(_, _, _, bodies) =>
          // add force contribution of each body by calling addForce
          bodies.foreach(b => addForce(b.mass, b.x, b.y))
        case a @ Fork(nw, ne, sw, se) => {
            val isSingle = a.size / distance(a.massX, a.massY, x, y) < theta
            if (isSingle) {
              addForce(a.mass, a.massX, a.massY)
            } else {
              traverse(nw)
              traverse(ne)
              traverse(sw)
              traverse(se)
            }
        }
          // see if node is far enough from the body,
          // or recursion is needed
      }

      traverse(quad)

      val nx = x + xspeed * delta
      val ny = y + yspeed * delta
      val nxspeed = xspeed + netforcex / mass * delta
      val nyspeed = yspeed + netforcey / mass * delta

      new Body(mass, nx, ny, nxspeed, nyspeed)
    }

  }

  val SECTOR_PRECISION = 8

  class SectorMatrix(val boundaries: Boundaries, val sectorPrecision: Int) {
    val sectorSize = boundaries.size / sectorPrecision
    val matrix = new Array[ConcBuffer[Body]](sectorPrecision * sectorPrecision)
    for (i <- 0 until matrix.length) matrix(i) = new ConcBuffer

    def +=(b: Body): SectorMatrix = {
      val x = (( maxMinX( b.x ) - boundaries.minX) / sectorSize ).toInt
      val y = (( maxMinY( b.y ) - boundaries.minY) / sectorSize ).toInt

      val i = sectorPrecision * y + x

      matrix(i) = matrix(i) += b

      this
    }

    def maxMinX(value: Float) = Math.min(Math.max(boundaries.minX, value), boundaries.maxX)
    def maxMinY(value: Float) = Math.min(Math.max(boundaries.minY, value), boundaries.maxY)

    def apply(x: Int, y: Int) = matrix(y * sectorPrecision + x)

//    def combine(that: SectorMatrix): SectorMatrix = {
//      val comb = matrix.zip(that.matrix).map { case (x, y) => x.combine(y) }
//
//      for (i <- 0 until matrix.length) matrix(i) = comb(i)
//
//      this
//    }

    def combine(that: SectorMatrix): SectorMatrix = {
      for (i <- 0 until matrix.length)
        this.matrix.update(i, this.matrix(i).combine(that.matrix(i)))
      this
    }

    def toQuad(parallelism: Int): Quad = {
      def BALANCING_FACTOR = 4
      def quad(x: Int, y: Int, span: Int, achievedParallelism: Int): Quad = {
        if (span == 1) {
          val sectorSize = boundaries.size / sectorPrecision
          val centerX = boundaries.minX + x * sectorSize + sectorSize / 2
          val centerY = boundaries.minY + y * sectorSize + sectorSize / 2
          val emptyQuad: Quad = Empty(centerX, centerY, sectorSize)
          val sectorBodies = this(x, y)
          sectorBodies.foldLeft(emptyQuad)(_ insert _)
        } else {
          val nspan = span / 2
          val nAchievedParallelism = achievedParallelism * 4
          val (nw, ne, sw, se) =
            if (parallelism > 1 && achievedParallelism < parallelism * BALANCING_FACTOR) parallel(
              quad(x, y, nspan, nAchievedParallelism),
              quad(x + nspan, y, nspan, nAchievedParallelism),
              quad(x, y + nspan, nspan, nAchievedParallelism),
              quad(x + nspan, y + nspan, nspan, nAchievedParallelism)
            ) else (
              quad(x, y, nspan, nAchievedParallelism),
              quad(x + nspan, y, nspan, nAchievedParallelism),
              quad(x, y + nspan, nspan, nAchievedParallelism),
              quad(x + nspan, y + nspan, nspan, nAchievedParallelism)
            )
          Fork(nw, ne, sw, se)
        }
      }

      quad(0, 0, sectorPrecision, 1)
    }

    override def toString = s"SectorMatrix(#bodies: ${matrix.map(_.size).sum})"
  }

  class TimeStatistics {
    private val timeMap = collection.mutable.Map[String, (Double, Int)]()

    def clear() = timeMap.clear()

    def timed[T](title: String)(body: =>T): T = {
      var res: T = null.asInstanceOf[T]
      val totalTime = /*measure*/ {
        val startTime = System.currentTimeMillis()
        res = body
        (System.currentTimeMillis() - startTime)
      }

      timeMap.get(title) match {
        case Some((total, num)) => timeMap(title) = (total + totalTime, num + 1)
        case None => timeMap(title) = (0.0, 0)
      }

      println(s"$title: ${totalTime} ms; avg: ${timeMap(title)._1 / timeMap(title)._2}")
      res
    }

    override def toString = {
      timeMap map {
        case (k, (total, num)) => k + ": " + (total / num * 100).toInt / 100.0 + " ms"
      } mkString("\n")
    }
  }
}
