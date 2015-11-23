#include "MSCKF.h"
#include "RK.h"
#include "MVG.h"

using namespace std;
using namespace Eigen;

MSCKF::MSCKF() {
    // 各种参数的初始化
    m_state_limit = 30;
    Matrix3d R_imu_to_cam;
    R_imu_to_cam << -Vector3d::UnitY(), -Vector3d::UnitX(), -Vector3d::UnitZ();
    m_q_imu_to_cam = HamiltonToJPL(Quaterniond(R_imu_to_cam));
    m_p_cam_in_imu << 0.0065, 0.0638, 0.0;

    m_cov_ng.setIdentity();
    m_cov_nwg.setIdentity();
    m_cov_na.setIdentity();
    m_cov_nwa.setIdentity();
}

void MSCKF::setNoiseCov(const Matrix3d &cov_ng, const Matrix3d &cov_nwg, const Matrix3d &cov_na, const Matrix3d &cov_nwa) {
    m_cov_ng = cov_ng;
    m_cov_nwg = cov_nwg;
    m_cov_na = cov_na;
    m_cov_nwa = cov_nwa;
}

void MSCKF::initialize(const JPL_Quaternion &q, const Vector3d &bg, const Vector3d &v, const Vector3d &ba, const Vector3d &p, double g) {
    m_q = q;
    m_bg = bg;
    m_v = v;
    m_ba = ba;
    m_p = p;
    m_g = g;
}

// [1] 中的 (9)(12)和(13)，bg 和 ba 是固定的，所以没有包含在这里，而PII和Phi都需要在相同的时间内进行积分所以包含在了这里
// StateVector := q, v, p, PII, Phi
typedef tuple<Vector4d, Vector3d, Vector3d, Matrix15d, Matrix15d> StateVector;

// StateVector 的数乘，用于 RK 进行积分
inline StateVector operator*(double s, const StateVector &v2) {
    return make_tuple(s*get<0>(v2), s*get<1>(v2), s*get<2>(v2), s*get<3>(v2), s*get<4>(v2));
}

// StateVector 的加法，用于 RK 进行积分
inline void operator+=(StateVector &v, const StateVector &v1) {
    get<0>(v) += get<0>(v1);
    get<1>(v) += get<1>(v1);
    get<2>(v) += get<2>(v1);
    get<3>(v) += get<3>(v1);
    get<4>(v) += get<4>(v1);
}

class MotionSystem {
public:
    MotionSystem(const MSCKF &filter, double t, const Vector3d &w, const Vector3d &a) :
        dt(t - filter.m_t_old),
        w0(filter.m_w_old - filter.m_bg),
        dw(w - filter.m_w_old),
        a0(filter.m_a_old - filter.m_ba),
        da(a - filter.m_a_old),
        filter(filter)
    {}

    StateVector operator() (double t, StateVector s) {
        StateVector r;
        double tt = t - filter.m_t_old;
        Vector3d wt = w0 + dw*tt / dt;
        Vector3d at = a0 + da*tt / dt;

        const JPL_Quaternion& sq = get<0>(s);
        const Vector3d& sv = get<1>(s);
        const Vector3d& sp = get<2>(s);
        const Matrix15d& sPII = get<3>(s);
        const Matrix15d& sPhi = get<4>(s);

        JPL_Quaternion& rq = get<0>(r);
        Vector3d& rv = get<1>(r);
        Vector3d& rp = get<2>(r);
        Matrix15d& rPII = get<3>(r);
        Matrix15d& rPhi = get<4>(r);

        Matrix3d cqt = JPL_CT(sq);
        Matrix3d neg_cross_w = -JPL_Cross(wt);
        Matrix3d cross_a = JPL_Cross(at);

        // (9)
        rq = 0.5*(JPL_Omega(wt)*sq);
        rv = cqt*at; rv.z() += filter.m_g;
        rp = sv;

        // (12)
        Matrix15d FxPII;
        multiplyFM(FxPII, neg_cross_w, -cqt, cross_a, sPII);
        rPII = FxPII + FxPII.transpose();
        rPII.block<3, 3>(0, 0) += filter.m_cov_ng;
        rPII.block<3, 3>(3, 3) += filter.m_cov_nwg;
        rPII.block<3, 3>(6, 6) += filter.m_cov_na;
        rPII.block<3, 3>(9, 9) += filter.m_cov_nwa;

        // (13)
        multiplyFM(rPhi, neg_cross_w, -cqt, cross_a, sPhi);

        return r;
    }
private:
    void multiplyFM(Matrix15d &Result, const Matrix3d &neg_cross_w, const Matrix3d &neg_cqt, const Matrix3d &cross_a, const Matrix15d &M) {
        Result.setZero();
        for (int i = 0; i < 15; i += 3) {
            Result.block<3, 3>(0, i) = neg_cross_w*M.block<3, 3>(0, i) - M.block<3, 3>(3, i);
            Result.block<3, 3>(6, i) = neg_cqt*(cross_a*M.block<3, 3>(0, i) + M.block<3, 3>(9, i));
            Result.block<3, 3>(12, i) = M.block<3, 3>(6, i);
        }
    }

    double dt;
    Vector3d w0, dw;
    Vector3d a0, da;
    const MSCKF &filter;
};

RungeKutta<double, StateVector> rk;

void MSCKF::propagate(double t, const Vector3d &w, const Vector3d &a) {
    static const Matrix15d eye = Matrix15d::Identity();
    Matrix15d Phi;

    if (m_has_old) {
        MotionSystem system(*this, t, w, a);
        tie(m_q, m_v, m_p, m_PII, Phi) = rk.integrate(system, tie(m_q, m_v, m_p, m_PII, eye), m_t_old, t);
        m_q = JPL_Normalize(m_q);
    }
    m_t_old = t;
    m_w_old = w;
    m_a_old = a;
    m_has_old = true;
}

void MSCKF::track(double t, const unordered_map<size_t, pair<size_t, Vector2d>> &matches) {
    //
    // 首先进行跟踪，找出当前丢失的 track
    //
    unordered_map<size_t, vector<Vector2d>> continued_tracks; // 其中包含了所有跟踪下来的 track
    size_t useful_state_length = 0; // 记录跟踪下来的 track 的最大长度
    for (const auto & f : matches) {
        if (m_tracks.count(f.second.first)>0) { // 这个特征跟踪到了已有的特征
            continued_tracks[f.first].swap(m_tracks.at(f.second.first)); // 将跟踪到的旧 track 放到 continued_tracks 中
        }
        useful_state_length = max(useful_state_length, continued_tracks[f.first].size());
        continued_tracks[f.first].push_back(f.second.second); // 将跟踪到的特征的位置加入到 track 中
    }

    vector<vector<Vector2d>> lost_tracks; // 其中包含了丢失的 track
    for (auto &t : m_tracks) {
        if (t.second.size() > 1) { // 当且仅当这个丢失的 track 包含了超过两个点，我们才考虑它，否则没有意义
            lost_tracks.emplace_back(move(t.second));
            useful_state_length = max(useful_state_length, lost_tracks.back().size());
        }
    }

    //
    // 去除老旧无用的 states
    //
    if (useful_state_length < m_states.size()) {
        size_t num_to_erase = m_states.size() - useful_state_length;
        m_states.erase(m_states.begin(), m_states.begin() + num_to_erase);
        m_PIC.erase(m_PIC.begin(), m_PIC.begin() + num_to_erase);
        m_PCC.erase(m_PCC.begin(), m_PCC.begin() + num_to_erase);
        for (auto &p : m_PCC) {
            p.erase(p.begin(), p.begin() + num_to_erase);
        }
    }

    //
    // 如果状态数依旧达到上限，在加入新的 state 前我们需要删除一些

    // 将删除过的 state 和 track 记录在这里
    vector<pair<Matrix3d, Vector3d>> new_states;
    unordered_map<size_t, vector<Vector2d>> new_tracks;
    vector<MatrixXd> new_PIC;
    vector<vector<MatrixXd>> new_PCC;

    // 将被删除的 track 记录在这里，它还要结合 m_state 进行 update
    unordered_map<size_t, vector<pair<Vector2d, size_t>>> jumping_tracks;
    if (m_states.size() == m_state_limit) {
        for (size_t i = 0; i < m_states.size(); ++i) {
            if (i % 3 != 1) {
                new_states.emplace_back(m_states[i]);
                new_PIC.push_back(m_PIC[i]);
                new_PCC.push_back(vector<MatrixXd>());
                for (size_t j = 0; j < m_PCC[i].size(); ++j) {
                    if (j % 3 != 1) {
                        new_PCC.back().push_back(m_PCC[i][j]);
                    }
                }
            }
        }
        for (auto& t : continued_tracks) { // 检查每个跟踪上的 track，从中去除相关联的特征组成独立的 track
            size_t tbegin = m_states.size() - (t.second.size() - 1);
            vector<Vector2d>& nt = new_tracks[t.first];
            vector<pair<Vector2d, size_t>>& jt = jumping_tracks[t.first];
            for (size_t i = tbegin; i < m_states.size(); ++i) {
                if (i % 3 != 1) {
                    nt.emplace_back(t.second[i - tbegin]);
                }
                else {
                    jt.emplace_back(t.second[i - tbegin], i);
                }
            }
            nt.emplace_back(t.second.back());
        }
    }

    //
    // 进行 EKF 更新
    // TODO: (17)~(33) here

    for (size_t i = 0; i < lost_tracks.size(); ++i) {
        Vector3d p = LinearLSTriangulation(lost_tracks[i], m_states);
    }
    for (auto &t : jumping_tracks) {
        if (t.second.size() > 1) {
            Vector3d p = LinearLSTriangulation(t.second, m_states);
        }
    }

    //
    // 完成状态扩充
    //

    if (m_states.size() == m_state_limit) {
        m_states.swap(new_states);
        m_tracks.swap(new_tracks);
        m_PIC.swap(new_PIC);
        m_PCC.swap(new_PCC);
    }
    else { // 如果我们没从中间删除任何状态，所有track都在continued_tracks中
        m_tracks.swap(continued_tracks);
    }

    // 计算当前相机姿态
    // (14)
    JPL_Quaternion q_world_to_cam = JPL_Multiply(m_q_imu_to_cam, m_q);
    Vector3d imu_to_cam_shift_in_world = JPL_CT(m_q)*m_p_cam_in_imu;
    Vector3d p_cam_in_world = m_p + imu_to_cam_shift_in_world;

    // (16)
    Matrix<double, 6, 15> Jc; 
    Jc.setZero();
    Jc.block<3, 3>(0, 0) = JPL_C(m_q_imu_to_cam);
    Jc.block<3, 3>(3, 0) = JPL_Cross(imu_to_cam_shift_in_world);
    Jc.block<3, 3>(3, 12).setIdentity();
    Matrix<double, 15, 6> JcT = Jc.transpose();
    
    // (15)
    m_PCC.push_back(vector<MatrixXd>());
    for (size_t i = 0; i < m_PIC.size(); ++i) {
        m_PCC.back().push_back((Jc*m_PIC[i]).transpose());
    }
    m_PCC.back().push_back(Jc*m_PII*JcT);
    m_PIC.push_back(m_PII*JcT);

    // 增加新的状态
    Matrix3d R_world_to_cam = JPL_C(q_world_to_cam);
    Vector3d T_world_to_cam = -R_world_to_cam*p_cam_in_world;
    m_states.emplace_back(R_world_to_cam, T_world_to_cam);
}
