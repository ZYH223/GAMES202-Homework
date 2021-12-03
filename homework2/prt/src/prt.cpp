#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/ray.h>
#include <filesystem/resolver.h>
#include <sh/spherical_harmonics.h>
#include <sh/default_image.h>
#include <Eigen/Core>
#include <fstream>
#include <random>
#include <stb_image.h>

NORI_NAMESPACE_BEGIN

namespace ProjEnv
{
    std::vector<std::unique_ptr<float[]>>
    LoadCubemapImages(const std::string &cubemapDir, int &width, int &height,
                      int &channel)
    {
        std::vector<std::string> cubemapNames{"negx.jpg", "posx.jpg", "posy.jpg",
                                              "negy.jpg", "posz.jpg", "negz.jpg"};
        std::vector<std::unique_ptr<float[]>> images(6);
        for (int i = 0; i < 6; i++)
        {
            std::string filename = cubemapDir + "/" + cubemapNames[i];
            int w, h, c;
            float *image = stbi_loadf(filename.c_str(), &w, &h, &c, 3);
            if (!image)
            {
                std::cout << "Failed to load image: " << filename << std::endl;
                exit(-1);
            }
            if (i == 0)
            {
                width = w;
                height = h;
                channel = c;
            }
            else if (w != width || h != height || c != channel)
            {
                std::cout << "Dismatch resolution for 6 images in cubemap" << std::endl;
                exit(-1);
            }
            images[i] = std::unique_ptr<float[]>(image);
            int index = (0 * 128 + 0) * channel;
            // std::cout << images[i][index + 0] << "\t" << images[i][index + 1] << "\t"
            //           << images[i][index + 2] << std::endl;
        }
        return images;
    }

    const Eigen::Vector3f cubemapFaceDirections[6][3] = {
        {{0, 0, 1}, {0, -1, 0}, {-1, 0, 0}},  // negx
        {{0, 0, 1}, {0, -1, 0}, {1, 0, 0}},   // posx
        {{1, 0, 0}, {0, 0, -1}, {0, -1, 0}},  // negy
        {{1, 0, 0}, {0, 0, 1}, {0, 1, 0}},    // posy
        {{-1, 0, 0}, {0, -1, 0}, {0, 0, -1}}, // negz
        {{1, 0, 0}, {0, -1, 0}, {0, 0, 1}},   // posz
    };

    float CalcPreArea(const float &x, const float &y)
    {
        return std::atan2(x * y, std::sqrt(x * x + y * y + 1.0));
    }

    float CalcArea(const float &u_, const float &v_, const int &width,
                   const int &height)
    {
        // transform from [0..res - 1] to [- (1 - 1 / res) .. (1 - 1 / res)]
        // ( 0.5 is for texel center addressing)
        float u = (2.0 * (u_ + 0.5) / width) - 1.0;
        float v = (2.0 * (v_ + 0.5) / height) - 1.0;

        // shift from a demi texel, mean 1.0 / size  with u and v in [-1..1]
        float invResolutionW = 1.0 / width;
        float invResolutionH = 1.0 / height;

        // u and v are the -1..1 texture coordinate on the current face.
        // get projected area for this texel
        float x0 = u - invResolutionW;
        float y0 = v - invResolutionH;
        float x1 = u + invResolutionW;
        float y1 = v + invResolutionH;
        float angle = CalcPreArea(x0, y0) - CalcPreArea(x0, y1) -
                      CalcPreArea(x1, y0) + CalcPreArea(x1, y1);

        return angle;
    }

    // template <typename T> T ProjectSH() {}

    template <size_t SHOrder>
    std::vector<Eigen::Array3f> PrecomputeCubemapSH(const std::vector<std::unique_ptr<float[]>> &images,
                                                    const int &width, const int &height,
                                                    const int &channel)
    {
        std::vector<Eigen::Vector3f> cubemapDirs;
        cubemapDirs.reserve(6 * width * height);
        for (int i = 0; i < 6; i++)
        {
            Eigen::Vector3f faceDirX = cubemapFaceDirections[i][0];
            Eigen::Vector3f faceDirY = cubemapFaceDirections[i][1];
            Eigen::Vector3f faceDirZ = cubemapFaceDirections[i][2];
            for (int y = 0; y < height; y++)
            {
                for (int x = 0; x < width; x++)
                {
                    float u = 2 * ((x + 0.5) / width) - 1;
                    float v = 2 * ((y + 0.5) / height) - 1;
                    Eigen::Vector3f dir = (faceDirX * u + faceDirY * v + faceDirZ).normalized();
                    cubemapDirs.push_back(dir);
                }
            }
        }
        constexpr int SHNum = (SHOrder + 1) * (SHOrder + 1);
        std::vector<Eigen::Array3f> SHCoeffiecents(SHNum);
        for (int i = 0; i < SHNum; i++)
            SHCoeffiecents[i] = Eigen::Array3f(0);
        float sumWeight = 0;
        for (int i = 0; i < 6; i++)
        {
            for (int y = 0; y < height; y++)
            {
                for (int x = 0; x < width; x++)
                {
                    // TODO: here you need to compute light sh of each face of cubemap of each pixel
                    // TODO: 此处你需要计算每个像素下cubemap某个面的球谐系数
                    Eigen::Vector3f dir = cubemapDirs[i * width * height + y * width + x];
                    int index = (y * width + x) * channel;
                    Eigen::Array3f Le(images[i][index + 0], images[i][index + 1],
                                      images[i][index + 2]);
                    auto deltaArea = CalcArea((float)x, (float)y, width, height);
                    //auto deltaArea = CalcArea(((float)x) / width, ((float)y) / height, width, height);
                    auto dirNormalized = dir.cast<double>().normalized();
                    // 将Lighting项投影到SH基函数上的过程可以通过黎曼积分来完成而避免了对球面空间采样
                    for (int l = 0; l <= SHOrder; l++)
                    {
                        for (int m = -l; m <= l; m++)
                        {
                            SHCoeffiecents[i] += Le * deltaArea * sh::EvalSH(l, m, dirNormalized);
                        }
                    }
                }
                std::cout << 100.0f * (float)(i * height + y) / (6 * height) << "%\r";
            }
        }
        return SHCoeffiecents;
    }
}

class PRTIntegrator : public Integrator
{
public:
    static constexpr int SHOrder = 2;
    static constexpr int SHCoeffLength = (SHOrder + 1) * (SHOrder + 1);

    enum class Type
    {
        Unshadowed = 0,
        Shadowed = 1,
        Interreflection = 2
    };

    PRTIntegrator(const PropertyList &props)
    {
        /* No parameters this time */
        m_SampleCount = props.getInteger("PRTSampleCount", 100);
        m_CubemapPath = props.getString("cubemap");
        auto type = props.getString("type", "unshadowed");
        if (type == "unshadowed")
        {
            m_Type = Type::Unshadowed;
        }
        else if (type == "shadowed")
        {
            m_Type = Type::Shadowed;
        }
        else if (type == "interreflection")
        {
            m_Type = Type::Interreflection;
            m_Bounce = props.getInteger("bounce", 1);
        }
        else
        {
            throw NoriException("Unsupported type: %s.", type);
        }
    }

    virtual void preprocess(const Scene *scene) override
    {

        // Here only compute one mesh
        const auto mesh = scene->getMeshes()[0];
        // Projection environment
        auto cubePath = getFileResolver()->resolve(m_CubemapPath);
        auto lightPath = cubePath / "light.txt";
        auto transPath = cubePath / "transport.txt";
        std::ofstream lightFout(lightPath.str());
        std::ofstream fout(transPath.str());
        int width, height, channel;
        std::vector<std::unique_ptr<float[]>> images =
            ProjEnv::LoadCubemapImages(cubePath.str(), width, height, channel);
        auto envCoeffs = ProjEnv::PrecomputeCubemapSH<SHOrder>(images, width, height, channel);
        m_LightCoeffs.resize(3, SHCoeffLength);
        for (int i = 0; i < envCoeffs.size(); i++)
        {
            lightFout << (envCoeffs)[i].x() << " " << (envCoeffs)[i].y() << " " << (envCoeffs)[i].z() << std::endl;
            m_LightCoeffs.col(i) = (envCoeffs)[i];
        }
        std::cout << "Computed light sh coeffs from: " << cubePath.str() << " to: " << lightPath.str() << std::endl;
        // Projection transport
        m_TransportSHCoeffs.resize(SHCoeffLength, mesh->getVertexCount());
        fout << mesh->getVertexCount() << std::endl;
        //m_Type = Type::Unshadowed;
        if (m_Type == Type::Unshadowed) std::cout << "Unshadowed" << std::endl;
        else if (m_Type == Type::Interreflection) std::cout << "Interreflection" << std::endl;
        else if (m_Type == Type::Shadowed) std::cout << "Shadowed" << std::endl;
        else  std::cout << "Unknown" << std::endl;
        std::cout << "Direct Illumination" << std::endl;
        for (int i = 0; i < mesh->getVertexCount(); i++)
        {
            const Point3f &v = mesh->getVertexPositions().col(i);
            const Normal3f &n = mesh->getVertexNormals().col(i);
            // shFunc的作用是定义Lighting Transport项，能够根据传入的theta和phi计算出传输项结果
            auto shFunc = [&](double phi, double theta) -> double {
                Eigen::Array3d d = sh::ToVector(phi, theta);
                const auto wi = Vector3f(d.x(), d.y(), d.z());
                const auto diffuse = /*(1/M_PI) */std::max(n.dot(wi), 0.0f);
                if (m_Type == Type::Unshadowed)
                {
                    // TODO: here you need to calculate unshadowed transport term of a given direction
                    // TODO: 此处你需要计算给定方向下的unshadowed传输项球谐函数值

                    // 仅考虑Diffuse的BRDF，直接返回漫反射结果，其传输项如下
                    return diffuse;
                }
                else
                {
                    // TODO: here you need to calculate shadowed transport term of a given direction
                    // TODO: 此处你需要计算给定方向下的shadowed传输项球谐函数值

                    // 考虑Shadowed Diffuse的情况，需要通过光线投射判断Visibility
                    Ray3f ray(v, wi.normalized());
                    if (scene->rayIntersect(ray)) return 0.0f;
                    else return diffuse;
                }
            };
            // 将传输项投影到SH基函数空间需要对球面空间进行采样
            /*int vindex[] { 5594, 5597, 5596, 5000, 5011, 5073 };
            int shindex[]{ 0, 2 };*/
            auto shCoeff = sh::ProjectFunction(SHOrder, shFunc, m_SampleCount);
            for (int j = 0; j < shCoeff->size(); j++)
            {
                m_TransportSHCoeffs.col(i).coeffRef(j) = (*shCoeff)[j];
                /*for (int tvi = 0; tvi < 6; tvi++)
                    for (int tsi = 0; tsi < 2; tsi++)
                        if (i == vindex[tvi] && j == shindex[tsi])
                            std::cout << "(v" << vindex[tvi] << ",sh" << shindex[tsi] << "):" << m_TransportSHCoeffs.col(i).coeffRef(j) << std::endl;*/
                //if (std::abs((*shCoeff)[j]) > 1.0f) std::cout << (*shCoeff)[j] << std::endl;
            }
            std::cout << 100.0f * (float)i / mesh->getVertexCount() << "%\r";
        }
        // TODO 目前还存在bug
        if (m_Type == Type::Interreflection)// 如果Interreflection启动，则计算间接光照
        {
            std::cout << "Indirect Illumination" << std::endl;
            // TODO: leave for bonus
            // 要实现Interreflection，就要通过两次以上的光线传输过程来计算光线传输项
            // 首先计算直接光照的球协系数，然后根据迭代求出间接光照的球协系数

            // 初始化迭代所需的传输项SH系数
            auto currentTransportSHCoeffs = Eigen::MatrixXf(m_TransportSHCoeffs);
            /*int vindex[]{ 5594, 5597, 5596, 5000, 5011, 5073 };
            int shindex[]{ 0, 2 };
            for (int tvi = 0; tvi < 6; tvi++)
                for (int tsi = 0; tsi < 2; tsi++)
                    std::cout << "(v" << vindex[tvi] << ",sh" << shindex[tsi] << "):" 
                    << m_TransportSHCoeffs.col(vindex[tvi]).coeffRef(shindex[tsi]) 
                    << "-" << currentTransportSHCoeffs.col(vindex[tvi]).coeffRef(shindex[tsi]) << std::endl;*/
            // 直接光照传输Td，n-bounce间接光照传输Ti(n)
            // 第n次迭代，用Td+Ti(n-1)生成Ti(n)
            const int maxBounces = 1;
            for (int bounce = 0; bounce < maxBounces; bounce++)
            {
                for (int i = 0; i < mesh->getVertexCount(); i++)
                {
                    // 跳过所有已经有直接光照结果的点，只计算直接光照中被遮挡的点
                    // 一般直接光照中visibility不为0的点是有SH系数结果的
                    if (std::abs(m_TransportSHCoeffs.col(i).coeffRef(0)) > 1e-8f) {
                        //std::cout << i << " skipped" << std::endl;
                        continue;
                    }
                    const Point3f& v = mesh->getVertexPositions().col(i);
                    const Normal3f& n = mesh->getVertexNormals().col(i);
                    // 定义用于计算要投影的原函数
                    auto shFunc = [&](double phi, double theta) -> double {
                        Eigen::Array3d d = sh::ToVector(phi, theta);
                        const auto wi = Vector3f(d.x(), d.y(), d.z());
                        Ray3f ray(v, wi);
                        Intersection its;
                        // 如果存着遮挡，则计算间接光照结果
                        if (scene->rayIntersect(ray, its))
                        {
                            float transportReflected = 0.0f;
                            const auto diffuse = std::max(n.dot(wi), 0.0f);
                            const auto bx = its.bary.x(); const auto by = its.bary.y(); const auto bz = its.bary.z();
                            const auto& tpx = currentTransportSHCoeffs.col(its.tri_index.x());
                            const auto& tpy = currentTransportSHCoeffs.col(its.tri_index.y());
                            const auto& tpz = currentTransportSHCoeffs.col(its.tri_index.z());
                            // 根据Td+Ti(n-1)的SH系数计算出被反射的传输项的结果
                            for (int l = 0; l <= SHOrder; l++) {
                                for (int m = -l; m <= l; m++) {
                                    auto shIndex = sh::GetIndex(l, m);
                                    auto shValue = sh::EvalSH(l, m, phi, theta);
                                    const auto coeff = bx * currentTransportSHCoeffs.col(its.tri_index.x()).coeffRef(shIndex)
                                        + by * currentTransportSHCoeffs.col(its.tri_index.y()).coeffRef(shIndex)
                                        + bz * currentTransportSHCoeffs.col(its.tri_index.z()).coeffRef(shIndex);
                                    transportReflected += coeff * shValue;

                                    /*if (std::abs(coeff) > 1.0f) {
                                        std::cout << "\tsh" << shIndex << ":" << coeff * shValue << "= ("
                                            << bx << "*" << tpx.coeffRef(shIndex) << "+" << by << "*" << tpy.coeffRef(shIndex) << "+" << bz << "*" << tpz.coeffRef(shIndex)
                                            << ")*" << shValue << std::endl;
                                        std::cout << "\t\t(v" << its.tri_index.x() << ",sh" << shIndex << "):" << tpx.coeffRef(shIndex)
                                            << "\t(v" << its.tri_index.y() << ",sh" << shIndex << "):" << tpy.coeffRef(shIndex)
                                            << "\t(v" << its.tri_index.z() << ",sh" << shIndex << "):" << tpz.coeffRef(shIndex) << std::endl;
                                    }*/
                                    
                                }
                            }
                            // 被反射的传输项，乘以cos项就是间接光照的传输项
                            /*if (std::abs(diffuse * transportReflected) > 2.0f)
                            std::cout << "indirect transport: " << diffuse * transportReflected << "=" << diffuse << "*" << transportReflected << std::endl;*/
                            return diffuse * transportReflected;
                        }
                        // 否则间接光照结果返回0（因为未遮挡说明没有间接光照）
                        // the light is direct then return 0, cause this is the computation for indirect shading
                        return 0;
                    };
                    // 将传输项投影到SH基函数空间需要对球面空间进行采样
                    auto shCoeff = sh::ProjectFunction(SHOrder, shFunc, m_SampleCount);
                    for (int j = 0; j < shCoeff->size(); j++)
                    {
                        // 此时的currentTransportSHCoeffs表示由Td+Ti(n-1)生成的Ti(n)
                        currentTransportSHCoeffs.col(i).coeffRef(j) = (*shCoeff)[j];
                    }
                    // 将Td+Ti(n)作为下一次迭代（生成Ti(n+1)）需要的SH系数
                    currentTransportSHCoeffs = m_TransportSHCoeffs + currentTransportSHCoeffs;

                    std::cout << 100.0f * (float)(i + bounce * mesh->getVertexCount()) / (maxBounces * mesh->getVertexCount()) << "%\r";
                }
            }
            m_TransportSHCoeffs = currentTransportSHCoeffs - m_TransportSHCoeffs;
        }

        // Save in face format
        for (int f = 0; f < mesh->getTriangleCount(); f++)
        {
            const MatrixXu &F = mesh->getIndices();
            uint32_t idx0 = F(0, f), idx1 = F(1, f), idx2 = F(2, f);
            for (int j = 0; j < SHCoeffLength; j++)
            {
                fout << m_TransportSHCoeffs.col(idx0).coeff(j) << " ";
            }
            fout << std::endl;
            for (int j = 0; j < SHCoeffLength; j++)
            {
                fout << m_TransportSHCoeffs.col(idx1).coeff(j) << " ";
            }
            fout << std::endl;
            for (int j = 0; j < SHCoeffLength; j++)
            {
                fout << m_TransportSHCoeffs.col(idx2).coeff(j) << " ";
            }
            fout << std::endl;
        }
        std::cout << "Computed SH coeffs"
                  << " to: " << transPath.str() << std::endl;
    }

    Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const
    {
        Intersection its;
        if (!scene->rayIntersect(ray, its))
            return Color3f(0.0f);

        const Eigen::Matrix<Vector3f::Scalar, SHCoeffLength, 1> sh0 = m_TransportSHCoeffs.col(its.tri_index.x()),
                                                                sh1 = m_TransportSHCoeffs.col(its.tri_index.y()),
                                                                sh2 = m_TransportSHCoeffs.col(its.tri_index.z());
        const Eigen::Matrix<Vector3f::Scalar, SHCoeffLength, 1> rL = m_LightCoeffs.row(0), gL = m_LightCoeffs.row(1), bL = m_LightCoeffs.row(2);
        
        auto dotValid = [&](Eigen::Matrix<Vector3f::Scalar, SHCoeffLength, 1> vec1, Eigen::Matrix<Vector3f::Scalar, SHCoeffLength, 1> vec2) {
            auto sum = 0.0f;
            for (auto i = 0; i < SHCoeffLength; i++)
            {
                sum += std::max(vec1[i] * vec2[i], 0.0f);
            }
            return sum;
        };

        Color3f c0 = Color3f(dotValid(rL, sh0), dotValid(gL, sh0), dotValid(bL, sh0)),
            c1 = Color3f(dotValid(rL, sh1), dotValid(gL, sh1), dotValid(bL, sh1)),
            c2 = Color3f(dotValid(rL, sh2), dotValid(gL, sh2), dotValid(bL, sh2));
        /*std::cout << "Color:"
            << std::endl << "c0:" << c0
            << std::endl << "c1:" << c1
            << std::endl << "c2:" << c2 << std::endl;*/
        //Color3f c0 = Color3f(rL.dot(sh0), gL.dot(sh0), bL.dot(sh0)),
        //        c1 = Color3f(rL.dot(sh1), gL.dot(sh1), bL.dot(sh1)),
        //        c2 = Color3f(rL.dot(sh2), gL.dot(sh2), bL.dot(sh2));

        const Vector3f &bary = its.bary;
        Color3f c = bary.x() * c0 + bary.y() * c1 + bary.z() * c2;
        /*std::cout << "Lighting:"
            << std::endl << "r:" << rL 
            << std::endl << "g:" << gL
            << std::endl << "b:" << bL << std::endl;
        std::cout << "Light Transport:"
            << std::endl << "sh0:" << sh0
            << std::endl << "sh1:" << sh1
            << std::endl << "sh2:" << sh2 << std::endl;*/
        /*std::cout << "Barycentric:" << bary << std::endl;
        std::cout << "Color:" << c << std::endl;*/
        // TODO: you need to delete the following four line codes after finishing your calculation to SH,
        //       we use it to visualize the normals of model for debug.
        // TODO: 在完成了球谐系数计算后，你需要删除下列四行，这四行代码的作用是用来可视化模型法线
        /*if (c.isZero()) {
            auto n_ = its.shFrame.n.cwiseAbs();
            return Color3f(n_.x(), n_.y(), n_.z());
        }*/
        return c;
    }

    std::string toString() const
    {
        return "PRTIntegrator[]";
    }

private:
    Type m_Type;
    int m_Bounce = 1;
    int m_SampleCount = 100;
    std::string m_CubemapPath;
    Eigen::MatrixXf m_TransportSHCoeffs;
    Eigen::MatrixXf m_LightCoeffs;
};

NORI_REGISTER_CLASS(PRTIntegrator, "prt");
NORI_NAMESPACE_END