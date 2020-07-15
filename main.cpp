#include "api.hpp"

Vec radiance(const Ray &r, int depth, unsigned short *Xi)
{
    double t;   // distance to intersection
    int id = 0; // id of intersected object
    if (!intersect(r, t, id))
        return Vec();                // if miss, return black
    const Sphere &obj = spheres[id]; // the hit object
    Vec x = r.o + r.d * t, n = (x - obj.p).norm(), nl = n.dot(r.d) < 0 ? n : n * -1, f = obj.c;
    double p = f.x > f.y && f.x > f.z ? f.x : f.y > f.z ? f.y : f.z; // max refl
    if (++depth > 5)
        if (erand48(Xi) < p)
            f = f * (1 / p);
        else
            return obj.e; //R.R.
    if (obj.refl == DIFF)
    { // Ideal DIFFUSE reflection
        double r1 = 2 * M_PI * erand48(Xi), r2 = erand48(Xi), r2s = sqrt(r2);
        Vec w = nl, u = ((fabs(w.x) > .1 ? Vec(0, 1) : Vec(1)) % w).norm(), v = w % u;
        Vec d = (u * cos(r1) * r2s + v * sin(r1) * r2s + w * sqrt(1 - r2)).norm();
        return obj.e + f.mult(radiance(Ray(x, d), depth, Xi));
    }
    else if (obj.refl == SPEC) // Ideal SPECULAR reflection
    {
        //f = f.mult(AiryTerm(650, std::abs(n.dot(r.d)), 1.0, 1.33, 1.0, 1.0, false));
        f = f.mult(AiryTerm(400, std::abs(n.dot(r.d)), 1.0, 1.7, 1.0, 0.0, false));
        return obj.e + f.mult(radiance(Ray(x, r.d - n * 2 * n.dot(r.d)), depth, Xi));
    }
    Ray reflRay(x, r.d - n * 2 * n.dot(r.d)); // Ideal dielectric REFRACTION
    Vec Rcolor = AiryTerm(400, std::abs(n.dot(r.d)), 1.0, 1.7, 1.0, 0.0, false);
    Vec Tcolor =  Vec(1,1,1) - Rcolor;
    double Pt = (Tcolor.x + Tcolor.y + Tcolor.z)/3.0;
    if(erand48(Xi) <= Pt)
    {
        f = f.mult(Tcolor) * (1/Pt);
        return obj.e + f.mult(radiance(Ray(x, r.d), depth, Xi));
    }
    else
    {
        f = f.mult(Rcolor) * (1/(1 - Pt));
        return obj.e + f.mult(radiance(reflRay, depth, Xi));   
    }
}
int main(int argc, char *argv[])
{
    int w = 1024, h = 768, samps = argc == 2 ? atoi(argv[1]) / 4 : 1; // # samples
    Ray cam(Vec(50, 52, 295.6), Vec(0, -0.042612, -1).norm());        // cam pos, dir
    Vec cx = Vec(w * .5135 / h), cy = (cx % cam.d).norm() * .5135, r, *c = new Vec[w * h];
#pragma omp parallel for schedule(dynamic, 1) private(r) // OpenMP
    for (int y = 0; y < h; y++)
    { // Loop over image rows
        fprintf(stderr, "\rRendering (%d spp) %5.2f%%", samps * 4, 100. * y / (h - 1));
        for (unsigned short x = 0, Xi[3] = {0, 0, y * y * y}; x < w; x++) // Loop cols
            for (int sy = 0, i = (h - y - 1) * w + x; sy < 2; sy++)       // 2x2 subpixel rows
                for (int sx = 0; sx < 2; sx++, r = Vec())
                { // 2x2 subpixel cols
                    for (int s = 0; s < samps; s++)
                    {
                        double r1 = 2 * erand48(Xi), dx = r1 < 1 ? sqrt(r1) - 1 : 1 - sqrt(2 - r1);
                        double r2 = 2 * erand48(Xi), dy = r2 < 1 ? sqrt(r2) - 1 : 1 - sqrt(2 - r2);
                        Vec d = cx * (((sx + .5 + dx) / 2 + x) / w - .5) +
                                cy * (((sy + .5 + dy) / 2 + y) / h - .5) + cam.d;
                        r = r + radiance(Ray(cam.o + d * 140, d.norm()), 0, Xi) * (1. / samps);
                    } // Camera rays are pushed ^^^^^ forward to start in interior
                    c[i] = c[i] + Vec(clamp(r.x), clamp(r.y), clamp(r.z)) * .25;
                }
    }
    FILE *f = fopen("image.ppm", "w"); // Write image to PPM file.
    fprintf(f, "P3\n%d %d\n%d\n", w, h, 255);
    for (int i = 0; i < w * h; i++)
        fprintf(f, "%d %d %d ", toInt(c[i].x), toInt(c[i].y), toInt(c[i].z));
}