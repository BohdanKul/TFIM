class Vertex
{
    public:
        Vertex(const float& _J12, const float& _h1, const float& _h2, const float _delta): J12(_J12), h1(_h1), h2(_h2), delta(_delta);
    private:
        float J12;
        float h1;
        float h2;
        float delta;

}


