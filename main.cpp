#include <algorithm>
#include <cfloat>
#include <cstdlib>
#include <fstream>
#include <iostream>
#define STB_IMAGE_IMPLEMENTATION
const double PI = 3.141592653;

#include "box.h"
#include "bvh.h"
#include "camera.h"
#include "hitablelist.h"
#include "jyorand.h"
#include "kuinkerm.h"
#include "material.h"
#include "perlin.h"
#include "rectangle.h"
#include "smoke.h"
#include "sphere.h"
#include "stb_image.h"
#include "onb.h"
#include "texture.h"
#include "transformation.h"

using namespace std;

Rand jyorandengine;
hitable_list world;

// 拒绝方法随机生成球内一点
vec3 randomInUnitSphere() {
  vec3 p;
  do {
    p = vec3(jyorandengine.jyoRandGetReal<double>(-1, 1),
             jyorandengine.jyoRandGetReal<double>(-1, 1),
             jyorandengine.jyoRandGetReal<double>(-1, 1));
  } while (p.squared_length() >= 1.0);
  return p;
}

// 拒绝方法随机生成半球内一点
vec3 randomInHemisphere(const vec3& normal) {
    // 随机生成一个单位球内的点
    vec3 in_unit_sphere = randomInUnitSphere(); 
    // 如果点在法线方向的半球上
    if (dot(in_unit_sphere, normal) > 0.0) 
        return in_unit_sphere; // 返回该点作为散射方向
    else
        return -in_unit_sphere; // 否则返回该点的反方向作为散射方向
}

// 拒绝方法随机生成圆内一点
vec3 randomInUnitDisk() {
  vec3 p;
  do {
    p = vec3(jyorandengine.jyoRandGetReal<double>(-1, 1),
             jyorandengine.jyoRandGetReal<double>(-1, 1), 0);
  } while (p.squared_length() >= 1.0);
  return p;
}

// 反演方法余弦采样生成半球面上一点
vec3 randomCosineDirection() {
    auto r1 = jyorandengine.jyoRandGetReal<double>(0,1);
    auto r2 = jyorandengine.jyoRandGetReal<double>(0,1);

    auto phi = 2*PI*r1;
    auto x = cos(phi)*sqrt(r2);
    auto y = sin(phi)*sqrt(r2);
    auto z = sqrt(1-r2);

    return vec3(x, y, z);
}

// 颜色着色
vec3 color(const ray& in, int depth) {
  hit_record rec;
  // 减少误差，-0.00001也可以是交点
  if (world.hitanythingbvh(in, 0.001, DBL_MAX, rec)) {
    // 反射出来的光线
    ray scattered;
    // 材料的吸收度
    vec3 attenuation;
    vec3 emitted = rec.mat_ptr->emitted(in, rec, rec.u, rec.v, rec.p);
    double pdf = -1.0;
    if (depth < 50 && rec.mat_ptr->scatter(in, rec, attenuation, scattered, pdf)) {
      // 强势硬编码，接管光线处理
      // 光源上随机一点
      vec3 on_light = vec3(jyorandengine.jyoRandGetReal<double>(213, 343), 554, jyorandengine.jyoRandGetReal<double>(227, 332));
      // 接触点到光源的射线
      vec3 to_light = on_light - rec.p;
      // 接触点到光源上一点的平方
      auto distance_squard = to_light.squared_length();
      // 单位化变成单位向量
      to_light = unit_vector(to_light);
      
      // 夹角为钝角，光源在物体背面，届不到
      if(dot(to_light, rec.normal) < 0)
        return emitted;

      double light_area = (343-213) * (332-227);
      // 单位向量巧妙计算夹角
      double light_cosine = fabs(to_light.y());
      // 如果光源和物体表面平行，也算届不到
      if(light_cosine < 1e-5)
        return emitted;

      pdf = distance_squard / (light_area * light_cosine);
      scattered = ray(rec.p, to_light, in.time());

      // 散射光线概率密度函数
      double scatteringpdf = rec.mat_ptr->scattering_pdf(in, rec, scattered);
      return emitted +  scatteringpdf * attenuation * color(scattered, depth + 1) / pdf;
    }
    else {
      // 直视光源则可以看到光源原本的颜色
      // if (!depth) emitted.make_unit_vector();
      return emitted;
    }
  } else {
    return vec3(0, 0, 0);
    // 天空
    // // 将射线in单位化，让其长度为1
    // vec3 unit_direction = unit_vector(in.direction());
    // // 单位化后的射线的y的取值范围为[-1,1]，将其+1再除以2转化为[0,1]
    // double t = 0.5 * (unit_direction.y() + 1.0);
    // // 进行颜色插值
    // return (1.0 - t) * vec3(1.0, 1.0, 1.0) + t * vec3(0.5, 0.7, 1.0);
  }
  exit(0);
}
std::vector<shared_ptr<hitable>> worldlist;
void buildWorld() {
  texture* whitelightptr = new constant_texture(vec3(17, 17, 17));
  texture* mikuptr = new constant_texture(vec3(0.223, 0.773, 0.733));
  texture* redptr = new constant_texture(vec3(0.65, 0.05, 0.05));
  texture* whiteptr = new constant_texture(vec3(0.73, 0.73, 0.73));
  texture* greenptr = new constant_texture(vec3(0.12, 0.45, 0.15));
  texture* groundtexptr = new constant_texture(vec3(0.48, 0.83, 0.53));

  texture* checkertextptr =
      new checker_texture(new constant_texture(vec3(0.2, 0.3, 0.1)),
                          new constant_texture(vec3(0.9, 0.9, 0.9)));

  texture* metaltextureptr = new constant_texture(
      vec3(0.5 * (1 + jyorandengine.jyoRandGetReal<double>(0, 1) *
                          jyorandengine.jyoRandGetReal<double>(0, 1)),
           0.5 * (1 + jyorandengine.jyoRandGetReal<double>(0, 1) *
                          jyorandengine.jyoRandGetReal<double>(0, 1)),
           0.5 * (1 + jyorandengine.jyoRandGetReal<double>(0, 1) *
                          jyorandengine.jyoRandGetReal<double>(0, 1))));
  texture* noisetextptr = new noise_texture(0.01);

  int nx, ny, nn;
  unsigned char* tex_data = stbi_load("earthmap.jpg", &nx, &ny, &nn, 0);
  texture* imagetextureptr = new image_texture(tex_data, nx, ny);

  worldlist.emplace_back(
      new rectangle_yz(0, 555, 0, 555, 555, new lambertian(redptr)));
  worldlist.emplace_back(
      new rectangle_yz(0, 555, 0, 555, 0, new lambertian(greenptr)));
  worldlist.emplace_back(new rectangle_xz(213, 343, 227, 332, 554,
                                          new diffuse_light(whitelightptr)));
  worldlist.emplace_back(
      new rectangle_xz(0, 555, 0, 555, 555, new lambertian(whiteptr)));
  worldlist.emplace_back(
      new rectangle_xz(0, 555, 0, 555, 0, new lambertian(whiteptr)));
  worldlist.emplace_back(
      new rectangle_xy(0, 555, 0, 555, 555, new lambertian(whiteptr)));
  worldlist.emplace_back(new translate(
      new rotate_y(
          new box(vec3(0, 0, 0), vec3(165, 165, 165), new lambertian(whiteptr)),
          -18),
      vec3(130, 0, 65)));
  worldlist.emplace_back(new translate(
      new rotate_y(
          new box(vec3(0, 0, 0), vec3(165, 330, 165), new lambertian(whiteptr)),
          15),
      vec3(265, 0, 295)));

  // 从世界列表中创建bvh树
  shared_ptr<hitable> rootptr;
  bvh_node(worldlist, rootptr);
  world = hitable_list(rootptr);
  // world = hitable_list(worldlist);
}

int getfileline() {
  std::ifstream file("output.PPM");
  // 判断文件是否打开成功
  if (file.is_open()) {
    int line_count = 0;
    std::string line;
    while (std::getline(file, line)) {
      ++line_count;
    }
    return line_count;
  } else
    return -1;
}

int main() {
  // 是否重新渲染
  int startoveragain = 1;

  int curline = getfileline();

  ofstream mout;
  if (startoveragain)
    mout.open("output.PPM");
  else
    mout.open("output.PPM", ios::app);

  // 画布的长
  int nx = 800;
  // 画布的宽
  int ny = 400;
  // 画布某一点的采样数量
  int ns = 10;

  buildWorld();
  vec3 lookfrom(278, 278, -800), lookat(278, 278, 0);
  camera cam(lookfrom, lookat, 40, double(nx) / double(ny), 0.0, 10.0, 0.0,
             1.0);

  int pauseflag = 1;
  int si, sj;
  if (curline <= 3 || startoveragain) {
    mout << "P3\n" << nx << " " << ny << "\n255\n";
    sj = ny - 1, si = 0;
  } else if (curline == nx * ny + 3) {
    system("pause");
    return 0;
  } else {
    curline -= 3;
    sj = ny - 1 - curline / nx;
    si = curline - nx * (curline / nx);
  }

  int sqrtns = int(sqrt(ns));
  double resqrtns = 1.0 / sqrtns;
  for (int j = sj; j >= 0; j--) {
    cout << "loading..." << 100 * (1.0 - double(j) / double(ny)) << "%";
    int starti = pauseflag ? si : 0;
    pauseflag = 0;
    for (int i = starti; i < nx; i++) {
      // 最终的颜色
      vec3 col(0, 0, 0);
      for (int dj = 0; dj < sqrtns; dj++)
        for (int di = 0; di < sqrtns; di++) {
          // 蒙特卡洛-抖动采样，将像素划分成更密的小格子，每个格子里随机取一个点采样
          double uplus =
              -0.5 + resqrtns * ((double)di +
                                 jyorandengine.jyoRandGetReal<double>(-1, 1));
          double vplus =
              -0.5 + resqrtns * ((double)dj +
                                 jyorandengine.jyoRandGetReal<double>(-1, 1));

          // 点(u,v)是点(i,j)的反离散化
          double u = (double(i) + uplus) / double(nx);
          double v = (double(j) + vplus) / double(ny);

          // 一条射向画布上点(u,v)的光线，注意(u,v)不是真实坐标而是在画布上的比例位置
          ray r = cam.get_ray(u, v);
          col += color(r, 0);
        }
      // 取颜色的平均值
      col /= double(ns);
      // gamma2修正，提升画面的质量
      col = vec3(sqrt(col[0]), sqrt(col[1]), sqrt(col[2]));
      int ir = int(255.99 * col[0]);
      int ig = int(255.99 * col[1]);
      int ib = int(255.99 * col[2]);
      ir = min(ir, 255), ig = min(ig, 255), ib = min(ib, 255);
      stringstream ss;
      ss << ir << " " << ig << " " << ib << "\n";
      mout << ss.str();
    }
    system("cls");
  }

  mout.close();
  system("pause");
  return 0;
}

// #include <cstdio>
// #include <cmath>
// #include <iostream>
// #include "vec3.h"
// #include "jyorand.h"

// using namespace std;

// Rand jyorandengine;
// const double PI = 3.141592653;

// inline vec3 random_cosine_direction() {
//     auto r1 = jyorandengine.jyoRandGetReal<double>(0,1);
//     auto r2 = jyorandengine.jyoRandGetReal<double>(0,1);

//     auto phi = 2*PI*r1;
//     auto x = cos(phi)*sqrt(r2);
//     auto y = sin(phi)*sqrt(r2);
//     auto z = sqrt(1-r2);

//     return vec3(x, y, z);
// }

// double f(const vec3& d) {
//     double cos_theta = d.z();
//     return cos_theta*cos_theta*cos_theta;
// }

// double q(const vec3& d) {
//     return d.z()/PI;
// }

// int main() {
//     int N = 10000000;
//     auto sum = 0.0;
//     for (int i = 0; i < N; i++) {
//         vec3 d = random_cosine_direction();
//         // 余弦分布的概率密度函数
//         sum += f(d) / q(d);
//     }
//     // 标准答案
//     std::cout << "PI/2 = " << PI / 2.0 << '\n';
//     // 均匀分布采样答案
//     std::cout << "Estimate = " << sum / N << '\n';
//     system("pause");
//     return 0;
// }